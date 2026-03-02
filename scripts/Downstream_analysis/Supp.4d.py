# -*- coding: utf-8 -*-
import os
import sys
from typing import Set

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('whitegrid')
    style.update(sns.axes_style('ticks'))
    sns.set(font_scale=1.4, style=style)

def load_gene_set(file_path: str) -> Set[str]:
    if not os.path.exists(file_path):
        print(f"Error: Gene set file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    with open(file_path, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def main():
    set_plotting_style()

    # --- Configuration ---
    cpm_matrix_path = "results/03_polya_matrix/combined_cpm_polya_matrix.csv"
    zygotic_genes_path = "data/annotation/zygotic_translation_genes_gonzalez_2025.txt"
    fig4_gene_sets_path = "results/downstream_analysis/12_fig4ab_scatters/fig4_gene_sets.csv"
    output_dir = "results/downstream_analysis/14_sperm_gene_sets"
    os.makedirs(output_dir, exist_ok=True)
    
    sample_to_analyze = 'sperm'
    min_count_num = 20 

    print("Loading data and gene sets...")
    merged_df = pd.read_csv(cpm_matrix_path)
    
    zygotic_translation_set = load_gene_set(zygotic_genes_path)
    
    fig4_sets_df = pd.read_csv(fig4_gene_sets_path)
    elongated_set = set(fig4_sets_df[fig4_sets_df['group'] == 'elongated']['gene_symbol'])
    consistent_set = set(fig4_sets_df[fig4_sets_df['group'] == 'consistent']['gene_symbol'])
    
    print(f"Filtering data for sample: {sample_to_analyze}...")
    read_count_col = f'{sample_to_analyze}_fulllength_polya_read_num'
    median_polya_col = f'{sample_to_analyze}_fulllength_median_polya'
    
    data_filtered = merged_df[merged_df[read_count_col] > min_count_num].copy()
    data_filtered.dropna(subset=[median_polya_col, 'gene_symbol'], inplace=True)

    def classify_gene_group(gene: str) -> str:
        if gene in zygotic_translation_set:
            return 'Zygotic Translation Set'
        elif gene in consistent_set:
            return 'Poly(A)-consistent'
        elif gene in elongated_set:
            return 'Poly(A)-elongated'
        else:
            return 'Other Genes'
    
    data_filtered['group'] = data_filtered['gene_symbol'].apply(classify_gene_group)
    
    print("Performing Mann-Whitney U tests...")
    zygotic_lengths = data_filtered[data_filtered['group'] == 'Zygotic Translation Set'][median_polya_col]
    consistent_lengths = data_filtered[data_filtered['group'] == 'Poly(A)-consistent'][median_polya_col]
    elongated_lengths = data_filtered[data_filtered['group'] == 'Poly(A)-elongated'][median_polya_col]

    p_vs_consistent = stats.mannwhitneyu(zygotic_lengths, consistent_lengths).pvalue
    p_vs_elongated = stats.mannwhitneyu(zygotic_lengths, elongated_lengths).pvalue
    print(f"Zygotic vs. Consistent P-value: {p_vs_consistent:.2e}")
    print(f"Zygotic vs. Elongated P-value: {p_vs_elongated:.2e}")

    print("Generating violin plot...")
    plt.figure(figsize=(8, 6))
    
    group_order = ['Poly(A)-consistent', 'Poly(A)-elongated', 'Zygotic Translation Set']
    palette = {'Poly(A)-consistent': 'skyblue', 'Poly(A)-elongated': 'salmon', 'Zygotic Translation Set': '#66c2a5'}

    sns.violinplot(
        x='group', y=median_polya_col, data=data_filtered,
        order=group_order, palette=palette, cut=0, bw_adjust=0.6
    )
    sns.stripplot(
        x='group', y=median_polya_col, data=data_filtered,
        order=group_order, color='black', size=2.5, jitter=0.2, alpha=0.4
    )

    y_max = data_filtered[median_polya_col].max()
    plt.plot([0, 0, 2, 2], [y_max*0.83, y_max*0.85, y_max*0.85, y_max*0.83], lw=1.5, c='black')
    plt.text(1, y_max*0.86, f'p = {p_vs_consistent:.2e}', ha='center', va='bottom', fontsize=12)
    
    plt.plot([1, 1, 2, 2], [y_max*0.78, y_max*0.8, y_max*0.8, y_max*0.78], lw=1.5, c='black')
    plt.text(1.5, y_max*0.81, f'p = {p_vs_elongated:.2e}', ha='center', va='bottom', fontsize=12)

    plt.title(f'Poly(A) Tail Length in {sample_to_analyze.capitalize()}', fontsize=16)
    plt.ylabel('Median Poly(A) Tail Length (nt)', fontsize=14)
    plt.xlabel('Gene Group', fontsize=14)
    plt.xticks(fontsize=12, rotation=0)
    plt.yticks(fontsize=12)
    
    output_path = os.path.join(output_dir, "sperm_gene_sets_polya_violin_plot.pdf")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    print(f"\nAnalysis complete. Plot saved to {output_path}")

if __name__ == "__main__":
    main()