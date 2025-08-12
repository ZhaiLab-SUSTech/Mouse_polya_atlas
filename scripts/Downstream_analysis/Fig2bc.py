# -*- coding: utf-8 -*-
import os
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import marsilea as ma
import marsilea.plotter as mp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))
from config_downstream import TISSUE_LIST, SHORT_TO_LONG_NAME

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('white')
    style.update(sns.axes_style('ticks'))
    sns.set(font_scale=1.4, style=style)

def plot_correlation_scatter(
    df: pd.DataFrame,
    x_sample: str,
    y_sample: str,
    min_count: int,
    output_dir: str
) -> None:
    query_str = f"`{x_sample}_fulllength_polya_read_num` > {min_count} and `{y_sample}_fulllength_polya_read_num` > {min_count}"
    merge_filter = df.query(query_str)
    
    if len(merge_filter) < 2:
        print(f"    - Skipping: Not enough common genes ({len(merge_filter)}) to plot.")
        return

    x_axis_col = f'{x_sample}_fulllength_median_polya'
    y_axis_col = f'{y_sample}_fulllength_median_polya'
    
    values = np.vstack([merge_filter[x_axis_col], merge_filter[y_axis_col]])
    try:
        kernel = stats.gaussian_kde(values)
        color_values = kernel(values)
    except np.linalg.LinAlgError:
        color_values = None

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(
        x=merge_filter[x_axis_col],
        y=merge_filter[y_axis_col],
        c=color_values,
        cmap="RdYlBu_r",
        s=1,
        vmin=0,
    )

    ax.plot([0, 350], [0, 350], color="black", ls="--", zorder=-10, linewidth=1)

    pearson, p_value = stats.pearsonr(merge_filter[x_axis_col], merge_filter[y_axis_col])
    ax.text(0.05, 0.95, f"r = {pearson:.2f}\nN = {len(merge_filter)}", 
            transform=ax.transAxes, fontsize=10, va='top', ha='left')

    ax.set_xlabel(f'{SHORT_TO_LONG_NAME.get(x_sample, x_sample)}\nMedian Poly(A) Length', fontsize=10)
    ax.set_ylabel(f'{SHORT_TO_LONG_NAME.get(y_sample, y_sample)}\nMedian Poly(A) Length', fontsize=10)
    sns.despine(ax=ax, top=True, right=True)

    ax.set_xlim(0, 350)
    ax.set_ylim(0, 350)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_yticks([0, 100, 200, 300])

    output_path = os.path.join(output_dir, f"scatter_{x_sample}_vs_{y_sample}.pdf")
    fig.savefig(output_path, bbox_inches='tight')
    plt.close(fig)

def main():
    set_plotting_style()

    # Configuration
    input_csv = "results/03_polya_matrix/combined_cpm_polya_matrix.csv"
    heatmap_output_dir = "results/downstream_analysis/03_correlation_heatmap"
    scatter_output_dir = "results/downstream_analysis/04_pairwise_scatters"
    min_read_count = 20
    
    pairs_to_plot: List[Tuple[str, str]] = [
        ('stomach', 'sperm'),
        ('testis', 'sperm'),
        ('adrenal', 'lung'),
        ('muscle', 'lung')
    ]
    
    if not os.path.exists(input_csv):
        print(f"Error: Input file not found at {input_csv}", file=sys.stderr)
        sys.exit(1)
    
    os.makedirs(heatmap_output_dir, exist_ok=True)
    os.makedirs(scatter_output_dir, exist_ok=True)
    merged_median_df = pd.read_csv(input_csv)

    read_num_cols = [col for col in merged_median_df.columns if col.endswith('_fulllength_polya_read_num')]
    mask = (merged_median_df[read_num_cols] > min_read_count).all(axis=1)
    merged_median_df_filtered = merged_median_df[mask]

    median_cols = [f"{tissue}_fulllength_median_polya" for tissue in TISSUE_LIST]
    median_cols_exist = [col for col in median_cols if col in merged_median_df_filtered.columns]
    
    merged_median_trim = merged_median_df_filtered[median_cols_exist]
    merged_median_trim.columns = [SHORT_TO_LONG_NAME.get(c.split('_')[0], c) for c in merged_median_trim.columns]

    corr_matrix = merged_median_trim.corr()

    h = ma.Heatmap(corr_matrix.values, cmap='Blues', width=6, height=6, vmax=1, vmin=0)
    plot_sample_names = mp.Labels(corr_matrix.columns.tolist(), fontsize=12, weight='bold')
    h.add_right(plot_sample_names, pad=0.1)
    h.add_bottom(plot_sample_names, pad=0.1)
    h.add_dendrogram('left', method="average", add_base=False, linewidth=1, colors='black', pad=0.1, size=1.3)
    h.add_dendrogram('top', method="average", add_base=False, linewidth=1, colors='black', pad=0.1, size=1.3)
    h.set_margin(.08)
    h.render()
    
    heatmap_output_path = os.path.join(heatmap_output_dir, "tissue_median_polya_correlation_heatmap.pdf")
    plt.savefig(heatmap_output_path, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to {heatmap_output_path}")

    for x_tissue, y_tissue in pairs_to_plot:
        plot_correlation_scatter(
            df=merged_median_df, 
            x_sample=x_tissue,
            y_sample=y_tissue,
            min_count=min_read_count,
            output_dir=scatter_output_dir
        )
    print(f"Scatter plots saved to {scatter_output_dir}")

if __name__ == "__main__":
    main()