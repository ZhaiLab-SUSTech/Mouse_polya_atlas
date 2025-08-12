# -*- coding: utf-8 -*-
import os
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))
from config_downstream import TISSUE_LIST, SHORT_TO_LONG_NAME

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('white')
    style.update(sns.axes_style('ticks'))
    style['xtick.major.size'] = 2
    style['ytick.major.size'] = 2
    sns.set(font_scale=1.4, style=style)

def generate_polya_bin(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = df.query('polya_length >= 10')
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()
        
    cut_bins = list(range(10, 350, 10))
    tmp_length_bin = np.digitize(df.polya_length, cut_bins)
    
    polya_bin = pd.crosstab(df.isoform_id, tmp_length_bin)
    polya_bin = polya_bin.divide(polya_bin.sum(axis=1), axis=0).fillna(0)
    
    mRNA_median = df.groupby('isoform_id').agg(
        count=('polya_length', 'size'),
        median=('polya_length', 'median')
    ).reset_index()
    
    bin_sort_value = polya_bin.apply(lambda x: np.argsort(-x.values), axis=1, result_type='expand')
    bin_sort_value.columns = [f"o{i}" for i in range(len(bin_sort_value.columns))]
    bin_sort_value['isoform_id'] = polya_bin.index
    
    bin_sort_value = bin_sort_value.merge(mRNA_median, on="isoform_id")
    bin_sort_value = bin_sort_value.sort_values(["o0", "median"])
    
    bin_sort_value = bin_sort_value.query('count >= 50')
    polya_bin = polya_bin.reindex(bin_sort_value.isoform_id).dropna()
    
    return polya_bin, bin_sort_value

def plot_polya_clustermap(df: pd.DataFrame, name: str, output_dir: str, vmax: float = 0.15) -> None:
    heatmap_width, heatmap_height, cbar_width = 3, 5, 0.2
    total_width, total_height = heatmap_width + cbar_width, heatmap_height

    cg = sns.clustermap(
        df, vmin=0, vmax=vmax,
        row_cluster=False, col_cluster=False,
        cmap="RdYlBu_r", figsize=(total_width, total_height),
        cbar_pos=None
    )
    
    cg.ax_heatmap.set_position([0.1, 0.1, heatmap_width / total_width, heatmap_height / total_height])
    cbar_pos = [0.1 + (heatmap_width / total_width) + 0.02, 0.1, 0.03, heatmap_height / total_height]
    cg.fig.add_axes(cbar_pos)
    plt.colorbar(cg.ax_heatmap.collections[0], cax=cg.fig.axes[-1])

    cg.ax_heatmap.tick_params(right=False, bottom=False, labelright=False, labelbottom=False)
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.set_ylabel("")

    plt.savefig(os.path.join(output_dir, f'{name.replace(" ", "_")}.png'), dpi=600, bbox_inches='tight')
    plt.close()

def main():
    set_plotting_style()
    
    # --- Configuration ---
    input_dir = "results/02_merged_readinfo"
    output_dir = "results/downstream_analysis/08_per_tissue_polya_heatmaps"
    os.makedirs(output_dir, exist_ok=True)
    
    polya_bins_list = []
    for tissue in tqdm(TISSUE_LIST, desc="Processing tissues"):
        try:
            rep1_path = os.path.join(input_dir, tissue, f"{tissue}_rep1.parquet")
            rep2_path = os.path.join(input_dir, tissue, f"{tissue}_rep2.parquet")
            
            s_rep1 = pd.read_parquet(rep1_path)
            s_rep2 = pd.read_parquet(rep2_path)
            
            polya_bin, _ = generate_polya_bin(pd.concat([s_rep1, s_rep2]))
            polya_bins_list.append(polya_bin)
        except FileNotFoundError:
            print(f"\nWarning: Could not find data for tissue '{tissue}'. Skipping.")
            polya_bins_list.append(pd.DataFrame())
            continue
            
    for i, tissue in enumerate(TISSUE_LIST):
        plot_name = SHORT_TO_LONG_NAME.get(tissue, tissue)
        polya_bin_df = polya_bins_list[i]
        
        if not polya_bin_df.empty:
            plot_polya_clustermap(polya_bin_df, name=plot_name, output_dir=output_dir)
        else:
            print(f"  - Skipping plot for {plot_name} (no data).")
            

if __name__ == "__main__":
    main()