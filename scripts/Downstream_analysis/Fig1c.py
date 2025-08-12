# -*- coding: utf-8 -*-
import os
import sys
from typing import List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import marsilea as ma
import marsilea.plotter as mp

# Import shared variables from the config file
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))
from config_downstream import TISSUE_LIST, SHORT_TO_LONG_NAME

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('white')
    sns.set(font_scale=1.4, style=style)

def create_distribution_df(data_series_list: List[pd.Series], min_len: int, max_len: int, bin_w: int) -> pd.DataFrame:
    bins = np.arange(min_len, max_len + bin_w, bin_w)
    dist_list = []
    
    for series in data_series_list:
        series_filtered = series[(series >= min_len) & (series < max_len)]
        total_count = len(series_filtered)
        if total_count == 0:
            dist_list.append(pd.Series(0, index=range(len(bins)-1)))
            continue

        digitized = np.digitize(series_filtered, bins)
        counts = pd.Series(digitized).value_counts()
        
        full_bins_df = pd.DataFrame({'bin': range(1, len(bins))})
        counts_df = counts.reset_index()
        counts_df.columns = ['bin', 'count']
        
        merged = full_bins_df.merge(counts_df, on='bin', how='left').fillna(0)
        
        normalized_counts = merged['count'] / total_count
        dist_list.append(normalized_counts)
        
    final_df = pd.concat(dist_list, axis=1)
    return final_df

def process_files_iteratively(
    base_dir: str, 
    tissues: List[str], 
    min_read_count: int
) -> Tuple[Dict[str, pd.Series], Dict[str, pd.Series]]:
    read_level_data = {}
    gene_median_data = {}

    for tissue in tissues:
        for i in [1, 2]:
            sample_name = f"{tissue}_rep{i}"
            file_path = os.path.join(base_dir, tissue, f"{sample_name}.parquet")
            
            if not os.path.exists(file_path):
                print(f"  - Warning: Could not find file {file_path}. Skipping.")
                continue

            print(f"  - Processing {sample_name}")
            df = pd.read_parquet(file_path)
            df.dropna(subset=['polya_length'], inplace=True)

            read_level_data[sample_name] = df['polya_length']

            gene_view = df.groupby('gene_id').agg(
                polya_length_median=('polya_length', 'median'),
                read_count=('polya_length', 'count')
            )
            filtered_medians = gene_view.query(f'read_count >= {min_read_count}')['polya_length_median']
            gene_median_data[sample_name] = filtered_medians

    print(f"Finished processing {len(read_level_data)} samples.")
    return read_level_data, gene_median_data

def main():
    set_plotting_style()
    
    input_dir = "results/02_merged_readinfo"
    output_dir = "results/downstream_analysis/02_summary_heatmaps"
    os.makedirs(output_dir, exist_ok=True)
    
    read_level_polya_data, gene_median_polya_data = process_files_iteratively(
        base_dir=input_dir,
        tissues=TISSUE_LIST,
        min_read_count=20
    )

    if not read_level_polya_data:
        print("No data found to process. Exiting.")
        return
        
    heatmap_order = [f"{t}_rep{i}" for t in TISSUE_LIST for i in [1, 2]]
    
    read_level_polya_list = [read_level_polya_data[sample] for sample in heatmap_order if sample in read_level_polya_data]
    gene_median_polya_list = [gene_median_polya_data[sample] for sample in heatmap_order if sample in gene_median_polya_data]
    
    # Plot Read-level Poly(A) Length Heatmap
    read_dist_df = create_distribution_df(read_level_polya_list, min_len=10, max_len=300, bin_w=1)
    read_dist_df.columns = [name for name in heatmap_order if name in read_level_polya_data]
    
    vsplit_idx = [i for i in range(2, len(read_dist_df.columns), 2)]
    
    h1 = ma.Heatmap(read_dist_df.T.values, cmap='Blues', width=4, height=12)
    h1.hsplit(cut=vsplit_idx)
    h1.add_left(mp.Chunk([SHORT_TO_LONG_NAME.get(t, t) for t in TISSUE_LIST], rotation=0, fontsize=12), pad=.05)
    h.add_dendrogram('left', colors=hex_colors, method="median",metric="euclidean",add_base=False, linewidth=0)
    h1.set_margin(.05)
    h1.render()
    plt.savefig(os.path.join(output_dir, "read_level_polya_heatmap.pdf"), bbox_inches='tight')
    plt.close()

    # Order From Dendrogram of Read-level Heatmap
    sorted_tissue_list = sorted_tissue_list = ['testis', 'sperm', 'pancreas', 'spleen', 'large', 'brain', 'thymus', 'thyroid', 'lung', 'kidney', 'adipose', 'heart', 'small', 'muscle', 'liver', 'stomach', 'adrenal', 'bone']
    # Plot Gene-median Poly(A) Length Heatmap
    print("Generating gene-median poly(A) length heatmap...")
    gene_dist_df = create_distribution_df(gene_median_polya_list, min_len=10, max_len=300, bin_w=5)
    gene_dist_df.columns = [name for name in heatmap_order if name in gene_median_polya_data]
    
    h2 = ma.Heatmap(gene_dist_df.T.values, cmap='Oranges', width=4, height=12)
    h2.hsplit(cut=vsplit_idx)
    h2.add_left(mp.Chunk([SHORT_TO_LONG_NAME.get(t, t) for t in sorted_tissue_list], rotation=0, fontsize=12), pad=.05)
    h2.set_margin(.05)
    h2.render()
    plt.savefig(os.path.join(output_dir, "gene_median_polya_heatmap.pdf"), bbox_inches='tight')
    plt.close()
    
    print(f"\nHeatmaps saved to {output_dir}")

if __name__ == '__main__':
    main()