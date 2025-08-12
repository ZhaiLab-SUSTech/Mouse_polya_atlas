# -*- coding: utf-8 -*-
import os
import sys
from typing import List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))
from config_downstream import TISSUE_LIST, SHORT_TO_LONG_NAME

def set_plotting_style():
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('white')
    style.update(sns.axes_style('ticks'))
    style['xtick.major.size'] = 2
    style['ytick.major.size'] = 2
    sns.set(font_scale=1.4, style=style)

def plot_overall_hist(
    samples: List[pd.DataFrame], 
    tissue_name: str, 
    output_path: str
) -> None:
    fig, ax = plt.subplots(figsize=(3, 1.2))
    
    min_length, max_length, bin_length = 10, 350, 1
    
    weights = [np.ones(len(s['polya_length'])) / len(s['polya_length']) for s in samples]
    labels = ['rep1', 'rep2']
    
    for i in range(2):
        plt.hist(
            samples[i]['polya_length'],
            bins=np.arange(min_length, max_length + bin_length, bin_length),
            weights=weights[i],
            histtype="step",
            label=labels[i],
            density=True, 
            linewidth=1.5,
            alpha=0.9
        )
        
    plt.ylim(0, 0.05)
    plt.xlim(min_length, max_length)
    plt.xticks([10, 50, 100, 150, 200, 250, 300, 350], fontsize=10)
    plt.yticks(fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    
    plt.title(SHORT_TO_LONG_NAME[tissue_name], y=0.92)
    plt.xlabel("Poly(A) tail length", fontsize=10)
    plt.ylabel("Ratio of Reads", fontsize=10)
    plt.legend(frameon=False, loc='upper right', prop={'size': 8})
    
    fig.savefig(os.path.join(output_path, f"{tissue_name}_polya_hist.pdf"), dpi=600, bbox_inches='tight')
    plt.close(fig)

def main():
    set_plotting_style()
    
    input_dir = "results/02_merged_readinfo"
    output_dir = "results/downstream_analysis/01_tissue_histograms"
    stats_file_path = os.path.join(output_dir, "tissue_read_stats.tsv")
    
    os.makedirs(output_dir, exist_ok=True)
    
    stats_records = []

    for tissue in TISSUE_LIST:
        try:
            rep1_path = os.path.join(input_dir, tissue, f"{tissue}_rep1.parquet")
            rep2_path = os.path.join(input_dir, tissue, f"{tissue}_rep2.parquet")
            
            samples = [
                pd.read_parquet(rep1_path).dropna(subset=['polya_length']),
                pd.read_parquet(rep2_path).dropna(subset=['polya_length'])
            ]
            
            plot_overall_hist(samples, tissue, output_dir)
            
            concat_sample = pd.concat(samples)
            gene_view = concat_sample.groupby('gene_id')['read_core_id'].nunique().reset_index(name='count')
            
            gt50_gene = (gene_view['count'] >= 50).sum()
            gt20_gene = (gene_view['count'] >= 20).sum()
            total_genes = len(gene_view)
            
            polya_num = (concat_sample['polya_length'] >= 10).sum()
            non_polya_num = (concat_sample['polya_length'] < 10).sum()
            
            stats_records.append({
                'tissue': tissue,
                'total_genes_detected': total_genes,
                'genes_with_over_20_reads': gt20_gene,
                'genes_with_over_50_reads': gt50_gene,
                'polya_reads_ge10': polya_num,
                'nonpolya_reads_lt10': non_polya_num
            })
            print(f"  - Processed {tissue}")

        except FileNotFoundError as e:
            print(f"Warning: Could not find data for tissue '{tissue}'. Skipping. Details: {e}")
            continue

    if stats_records:
        stats_df = pd.DataFrame(stats_records)
        stats_df.to_csv(stats_file_path, sep='\t', index=False)
        print(f"\nStatistics summary saved to {stats_file_path}")

if __name__ == '__main__':
    main()