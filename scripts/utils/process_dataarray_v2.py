#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Description: return raw polya list instead of bins as data based on process_dataarray.py
'''

import argparse
from collections import Counter, defaultdict
from multiprocessing import Pool
import pickle
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np
import scipy.stats
from sklearn.preprocessing import normalize
import os
import sys
import logging
from merge_utils import parse_gtf, merge_dict

# Configure logging to write to both console and a log file
def configure_logging(log_file_path: str):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path),
            logging.StreamHandler(sys.stdout)
        ]
    )

def load_concat_samples(merged_data_paths: List[str], min_polya: int, max_polya: int) -> Tuple[pd.DataFrame, List[str]]:
    samples = []
    samples_names = []
    for path in merged_data_paths:
        name = os.path.basename(path).split('.')[0] # get name like sample_rep1
        samples_names.append(name)
        logging.info(f'Start processing {name}')
        samples.append(parse_sample(path, name, min_polya, max_polya))
    return (pd.concat(samples), samples_names)

def parse_sample(path: str, name: str, min_polya: int, max_polya: int) -> pd.DataFrame:
    if not os.path.exists(path):
        logging.info(f"Error: Merged data file {path} does not exist.")
        sys.exit(1)
    try:
        merged = pd.read_parquet(path)
        merged.loc[:, 'sample_info'] = name
        logging.info(f"Loaded merged data from {path}")
    except Exception as e:
        logging.info(f"Error loading merged data from {path}: {e}")
        sys.exit(1)
    return merged

def generate_bins(df: pd.DataFrame, bin_num: np.ndarray):
    hist, _ = np.histogram(df.polya_length, bins=bin_num)
    normalized_hist = normalize(hist.reshape(1, -1), norm='l1')
    return normalized_hist
    

def process_gene(df: pd.DataFrame, gene: str, bin_num: np.ndarray, tissue: str):
    try:
        # return compute_statisitcs(df.query('new_isoform == @gene'))
        return df.query('new_isoform == @gene').polya_length.tolist()
    except Exception as e:
        logging.info(f"Error in processing gene {gene} in {tissue}: {e}")
        return None

def worker(df: pd.DataFrame, min_count: int, tissue: str, bin_num: np.ndarray, thread: int = 2) -> Tuple[List[str], List[int], List[np.ndarray], List[str]]:
    df['new_isoform'] = np.where(df['isoform_id']!='nan', df['isoform_id'], np.where(df['mRNA']!='nan', df['mRNA'], np.nan))
    gene_view = df.loc[:, ['new_isoform', 'polya_length']].groupby('new_isoform').agg(['count', 'median'])
    gene_view.columns = ['count', 'median']
    gene_list = gene_view.query('count>@min_count').index.tolist()

    tissue_info = [tissue] * len(gene_list)
    expression_info = gene_view.query('count>@min_count')['count'].tolist()

    logging.info(f'Processing {tissue}')

    with Pool(processes=thread) as pool:
        data_array = list(pool.starmap(process_gene, [(df, gene, bin_num, tissue) for gene in gene_list]))

    logging.info('Multiprocessing finished')
    # Remove None values from data_array
    # data_array = [data for data in data_array if data is not None]

    return (tissue_info, expression_info, data_array, gene_list)



def main():
    parser = argparse.ArgumentParser(description="Analyze merged sample data.")
    parser.add_argument('--merged_data', nargs='+', required=True, help='Path to the merged Parquet file(s)')
    parser.add_argument('--label', required=True, help='label for samples')
    parser.add_argument('--output_dir', required=True, help='Directory to save analysis results')
    parser.add_argument('--gtf_path', nargs='+', required=True, help='Path to gtf file')
    parser.add_argument('--min_count', type=int, default=20, help='Min read count for one gene')
    parser.add_argument('--bin_width', type=int, default=5, help='Histogram bin width')
    parser.add_argument('--thread', type=int, default=2, help='Thread for parallel running')

    args = parser.parse_args()

    merged_data_paths = args.merged_data
    output_dir = args.output_dir
    label = args.label
    gtf_paths = args.gtf_path
    min_count = args.min_count
    bin_width = args.bin_width
    thread = args.thread

    # histogram parameters
    min_num = 10
    max_num = 300

    # Configure logging
    log_file_path = os.path.join(output_dir, f'analyze_isoform_{label}.log')
    configure_logging(log_file_path)

    logging.info("Loading samples ...")
    rep_merged_df, samples_names = load_concat_samples(merged_data_paths, min_polya=min_num, max_polya=max_num )

    # parse gtf
    logging.info(f'Processing gtf files ...')
    _, transcript_to_gene1, gene_to_symbol1 = parse_gtf(gtf_paths[0])
    _, transcript_to_gene2, gene_to_symbol2 = parse_gtf(gtf_paths[1])

    transcript_to_gene = merge_dict(transcript_to_gene1, transcript_to_gene2)
    gene_to_symbol = merge_dict(gene_to_symbol1, gene_to_symbol2)

    bin_num = np.arange(min_num, max_num + bin_width, bin_width)

    tissue_info, expression_info, data_array, gene_list = worker(rep_merged_df, min_count=min_count, tissue=label, bin_num=bin_num, thread=int(thread))

    with open(f'{output_dir}/{label}_dataarray.pkl', 'wb') as file:
        pickle.dump(data_array, file)
    with open(f'{output_dir}/{label}_sampleinfo.pkl', 'wb') as file:
        pickle.dump(tissue_info, file)
    with open(f'{output_dir}/{label}_geneinfo_isoform.pkl', 'wb') as file:
        pickle.dump(gene_list, file)
    with open(f'{output_dir}/{label}_expression.pkl', 'wb') as file:
        pickle.dump(expression_info, file)


    gene_list_geneid = [transcript_to_gene[gene] if gene in transcript_to_gene.keys() else np.nan for gene in gene_list]
    with open(f'{output_dir}/{label}_geneinfo_geneid.pkl', 'wb') as file:
        pickle.dump(gene_list_geneid, file)

    gene_list_genesymbol = [gene_to_symbol[gene] if gene in gene_to_symbol.keys() else np.nan for gene in gene_list_geneid]
    with open(f'{output_dir}/{label}_geneinfo_genesymbol.pkl', 'wb') as file:
        pickle.dump(gene_list_genesymbol, file)

if __name__ == "__main__":
    main()




