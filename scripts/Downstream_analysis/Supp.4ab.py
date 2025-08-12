# -*- coding: utf-8 -*-
import os
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('white')
    style.update(sns.axes_style('ticks'))
    style['xtick.major.size'] = 2
    style['ytick.major.size'] = 2
    sns.set(font_scale=1.4, style=style)

def plot_expression_scatter(
    df: pd.DataFrame, 
    x_sample: str, 
    y_sample: str, 
    title: str, 
    output_path: str
) -> None:
    df_copy = df.copy()
    x_cpm_col = f'{x_sample}_bam_cpm'
    y_cpm_col = f'{y_sample}_bam_cpm'
    x_log_col = f'log_{x_cpm_col}'
    y_log_col = f'log_{y_cpm_col}'

    df_copy[x_log_col] = np.log2(df_copy[x_cpm_col] + 1) # Add 1 to avoid log(0)
    df_copy[y_log_col] = np.log2(df_copy[y_cpm_col] + 1)

    values = np.vstack([df_copy[x_log_col], df_copy[y_log_col]])
    kernel = stats.gaussian_kde(values)
    color_values = kernel(values)

    fig, ax = plt.subplots(figsize=(3, 3))
    ax.scatter(
        x=df_copy[x_log_col], y=df_copy[y_log_col], c=color_values,
        cmap="RdYlBu_r", s=8
    )
    ax.plot([0, 12], [0, 12], color="black", ls="--", zorder=-10)
    ax.set_xlabel(f'{x_sample}\nlog2(CPM+1)', fontsize=12)
    ax.set_ylabel(f'{y_sample}\nlog2(CPM+1)', fontsize=12)
    sns.despine(ax=ax, top=True, right=True)
    ax.set_xlim(2, 12)
    ax.set_ylim(2, 12)
    ax.set_xticks([2, 4, 6, 8, 10, 12])
    ax.set_yticks([2, 4, 6, 8, 10, 12])
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    set_plotting_style()
    
    # --- Configuration ---
    cpm_matrix_path = "results/03_polya_matrix/combined_cpm_polya_matrix.csv"
    output_dir = "results/downstream_analysis/12_fig4ab_scatters"
    os.makedirs(output_dir, exist_ok=True)

    gene_sets_output_path = os.path.join(output_dir, "fig4_gene_sets.csv")
    
    min_read_count = 20
    x_sample = 'stomach'
    y_sample = 'sperm'
    
    print(f"Loading data from {cpm_matrix_path}...")
    merged_df = pd.read_csv(cpm_matrix_path)
    
    query_str = f"`{x_sample}_fulllength_polya_read_num` > {min_read_count} and `{y_sample}_fulllength_polya_read_num` > {min_read_count}"
    merge_filter = merged_df.query(query_str).copy()

    x_polya_col = f'{x_sample}_fulllength_median_polya'
    y_polya_col = f'{y_sample}_fulllength_median_polya'
    
    values = np.vstack([merge_filter[x_polya_col], merge_filter[y_polya_col]])
    kernel = stats.gaussian_kde(values)
    color_values = kernel(values)
    
    high_density_points = merge_filter[color_values > 0.00015].copy()
    high_density_colors = color_values[color_values > 0.00015]

    fig, ax = plt.subplots(figsize=(4, 4))
    scatter = ax.scatter(
        x=high_density_points[x_polya_col], y=high_density_points[y_polya_col],
        c=high_density_colors, cmap="RdYlBu_r", s=1, vmin=0, vmax=0.0003
    )
    
    ax.plot([0, 350], [0, 350], color="black", ls="--", zorder=-10)
    ax.set_xlabel(f'{x_sample}\nMedian Poly(A) Length', fontsize=12)
    ax.set_ylabel(f'{y_sample}\nMedian Poly(A) Length', fontsize=12)
    sns.despine(ax=ax, top=True, right=True)
    ax.set_xlim(0, 350); ax.set_ylim(0, 350)
    ax.set_xticks([0, 100, 200, 300]); ax.set_yticks([0, 100, 200, 300])
    
    ax.hlines(y=135, xmin=80, xmax=150, color='r')
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label('Density', fontsize=15)
    
    fig.savefig(os.path.join(output_dir, "fig4a_polya_scatter.pdf"), bbox_inches='tight')
    plt.close()

    elongated_genes_df = high_density_points.query(f'`{y_polya_col}` > 135')
    consistent_genes_df = high_density_points.query(f'`{y_polya_col}` <= 135')
    print(f"Defined gene sets: {len(elongated_genes_df)} elongated, {len(consistent_genes_df)} consistent.")

    elongated_genes_df['group'] = 'elongated'
    consistent_genes_df['group'] = 'consistent'
    
    gene_sets_to_save = pd.concat([
        elongated_genes_df[['gene_symbol', 'group']],
        consistent_genes_df[['gene_symbol', 'group']]
    ])
    
    # Save to CSV
    gene_sets_to_save.to_csv(gene_sets_output_path, index=False)
    
    # ===========================================================
    plot_expression_scatter(
        df=elongated_genes_df, x_sample=x_sample, y_sample=y_sample,
        title="Elongated Genes",
        output_path=os.path.join(output_dir, "fig4b_expr_scatter_elongated.pdf")
    )
    plot_expression_scatter(
        df=consistent_genes_df, x_sample=x_sample, y_sample=y_sample,
        title="Consistent Genes",
        output_path=os.path.join(output_dir, "fig4b_expr_scatter_consistent.pdf")
    )
    
    print(f"\nAnalysis complete. Plots for Fig 4a and 4b saved to {output_dir}")

if __name__ == "__main__":
    main()