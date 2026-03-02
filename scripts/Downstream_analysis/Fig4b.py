# -*- coding: utf-8 -*-
import os
import sys
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))
from config_downstream import TISSUE_LIST

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('whitegrid')
    style.update(sns.axes_style('ticks'))
    sns.set(font_scale=1.2, style=style)

def calculate_gene_tissue_correlation(
    isoform_id: str, 
    df: pd.DataFrame,
    tissues: List[str]
) -> Tuple[Optional[float], Optional[float]]:
    try:
        gene_data = df.loc[isoform_id]
    except KeyError:
        return None, None

    median_polya_values = []
    cpm_values = []

    for tissue in tissues:
        polya_col = f'{tissue}_fulllength_median_polya'
        cpm_col = f'{tissue}_bam_cpm'
        
        if polya_col in gene_data and cpm_col in gene_data and \
           pd.notna(gene_data[polya_col]) and pd.notna(gene_data[cpm_col]):
            median_polya_values.append(gene_data[polya_col])
            cpm_values.append(gene_data[cpm_col])

    if len(median_polya_values) < 3:
        return None, None

    try:
        corr, p_value = stats.spearmanr(median_polya_values, cpm_values)
        return corr, p_value
    except ValueError:
        return None, None

def main():
    set_plotting_style()

    # --- Configuration ---
    cpm_matrix_path = "results/03_polya_matrix/combined_cpm_polya_matrix.csv"
    
    output_dir = "results/downstream_analysis/11_polya_expression_correlation"
    os.makedirs(output_dir, exist_ok=True)

    cpm_matrix = pd.read_csv(cpm_matrix_path)
    cpm_matrix.set_index('isoform_id', inplace=True)
    
    high_expression_isoforms_path = "results/downstream_analysis/high_expression_isoforms.txt"
    if not os.path.exists(high_expression_isoforms_path):
        print(f"Error: High expression isoform list not found.", file=sys.stderr)
        print(f"Please run 'extract_high_expression_isoforms.py' first.", file=sys.stderr)
        sys.exit(1)
    with open(high_expression_isoforms_path, 'r') as f:
        high_expression_genes = [line.strip() for line in f]

    corr_results = []
    for isoform in tqdm(high_expression_genes, desc="Processing isoforms"):
        corr, p = calculate_gene_tissue_correlation(isoform, cpm_matrix, TISSUE_LIST)
        if corr is not None:
            corr_results.append({'isoform_id': isoform, 'corr': corr, 'p_value': p})

    corr_abund_df = pd.DataFrame(corr_results)
    corr_abund_df['significance'] = ['p < 0.05' if p < 0.05 else 'p >= 0.05' for p in corr_abund_df['p_value']]
    
    output_csv_path = os.path.join(output_dir, "polya_cpm_correlation_results.csv")
    corr_abund_df.to_csv(output_csv_path, index=False)
    print(f"Correlation results saved to {output_csv_path}")

    plt.figure(figsize=(5, 4))
    sns.histplot(
        data=corr_abund_df, 
        x='corr', 
        hue='significance',
        multiple="stack", 
        binwidth=0.1,
        palette={'p < 0.05': '#d62728', 'p >= 0.05': 'darkgrey'},
        edgecolor="white",
        alpha=0.8
    )
    plt.title("Poly(A) Length vs. Expression Correlation", fontsize=14)
    plt.xlabel("Spearman Correlation (r)")
    plt.ylabel("Number of Isoforms")
    plt.xlim(-1, 1)
    
    output_plot_path = os.path.join(output_dir, "polya_cpm_correlation_histogram.pdf")
    plt.savefig(output_plot_path, bbox_inches='tight')
    plt.close()
    
    print(f"Final histogram saved to {output_plot_path}")

if __name__ == "__main__":
    main()