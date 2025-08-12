#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from collections import defaultdict
from typing import Dict, List, Tuple
import pandas as pd
import os
import sys
import logging


# Function to generate transcript_id to exon count mapping
def generate_transcript_exon_count(gtf_file: str) -> Dict[str, int]:
    exon_count_map = defaultdict(int)
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue  # Skip lines that are not exons

            attributes = fields[8]
            attribute_dict = {}
            for attr in attributes.split(';'):
                if attr.strip():
                    key_value = attr.strip().split(' ')
                    if len(key_value) >= 2:
                        key = key_value[0]
                        value = ' '.join(key_value[1:]).strip('"')
                        attribute_dict[key] = value

            transcript_id = attribute_dict.get('transcript_id', '')
            if transcript_id:
                exon_count_map[transcript_id] += 1

    return dict(exon_count_map)

def parse_sample(path: str, name: str, gtf_path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        logging.info(f"Error: Merged data file {path} does not lesexist.")
        sys.exit(1)
    try:
        merged = pd.read_parquet(path)
        merged = merged.query("assignment_type.str.contains('unique') and not assignment_events.str.contains('ism')")
        merged.loc[:, 'sample_info'] = name
        transcript_exon_count_map = generate_transcript_exon_count(gtf_path)
        merged.loc[:, 'exon_num'] = merged['isoform_id'].map(transcript_exon_count_map)
        merged.loc[:, 'exon_num'] = merged['exon_num'].astype(int)
        merged = merged.query('not ( (exon_num > 1) and assignment_events.str.contains("mono") )')
        logging.info(f"Loaded merged data from {path}")
    except Exception as e:
        logging.info(f"Error loading merged data from {path}: {e}")
        sys.exit(1)
    return merged

def process_gene_info(merge_df: pd.DataFrame) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    def generate_gene_isoform_mapping(df):
        gene_isoform_map = df.groupby('gene_id')['isoform_id'].unique().apply(list).to_dict()
        return gene_isoform_map

    def generate_gene_symbol_mapping(df):
        gene_symbol_map = df.drop_duplicates('gene_id').set_index('gene_id')['gene_symbol'].to_dict()
        return gene_symbol_map
    
    return generate_gene_isoform_mapping(merge_df), generate_gene_symbol_mapping(merge_df)