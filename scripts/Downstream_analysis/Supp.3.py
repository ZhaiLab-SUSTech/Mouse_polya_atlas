# -*- coding: utf-8 -*-
import os
import sys
from typing import List

import numpy as np
import pandas as pd
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

def plot_stacked_contribution(df: pd.DataFrame, name: str, output_dir: str) -> None:
    top_genes = df.query('isoform_id.str.startswith("EN")').groupby("gene_symbol").size().nlargest(10).index
    df["group"] = df["gene_symbol"].apply(lambda x: x if x in top_genes else "Other")
    
    grouped_data = [df[df["group"] == gene]["polya_length"] for gene in top_genes]
    grouped_data.append(df[df["group"] == "Other"]["polya_length"])
    
    colors = sns.color_palette("tab10", 10) + ["gray"]
    bin_edges = np.arange(10, 300 + 5, 5)

    plt.figure(figsize=(8, 6))
    plt.hist(
        grouped_data, bins=bin_edges, stacked=True, color=colors, density=True, 
        label=list(top_genes) + ["Other"], alpha=0.8, linewidth=0.5
    )
    plt.hist(
        df["polya_length"], bins=bin_edges, color="black", 
        histtype="step", linewidth=1.5, label="Total", density=True
    )

    plt.xlabel("PolyA Length")
    plt.ylabel("Ratio")
    plt.xlim(0, 300)
    plt.ylim(0, 0.015)
    plt.title(name)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    
    output_path = os.path.join(output_dir, f"contribution_stacked_hist_{name}.svg")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def plot_bar_genecount(df: pd.DataFrame, name: str, output_dir: str) -> None:
    gene_counts = df.query('isoform_id.str.startswith("EN")').groupby("gene_symbol").size().nlargest(10)
    colors = sns.color_palette("tab10", 10)[::-1]

    plt.figure(figsize=(3, 6))
    gene_counts.sort_values().plot.barh(
        color=colors,
        width=0.8
    )

    plt.xlabel("Transcript Count", fontsize=12)
    plt.ylabel(None)
    plt.xticks(fontsize=8)
    plt.gca().invert_yaxis()
    plt.xlim(0, gene_counts.max() * 1.15)
    
    output_path = os.path.join(output_dir, f"contribution_bar_chart_{name}.svg")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def main():
    set_plotting_style()

    # --- Configuration ---
    samples_to_process: List[str] = ['muscle_rep1', 'muscle_rep2', 'spleen_rep1', 'spleen_rep2']
    
    input_dir_base = "results/02_merged_readinfo"
    output_dir = "results/downstream_analysis/16_gene_contribution"
    os.makedirs(output_dir, exist_ok=True)
    
    for sample_name in samples_to_process:
        tissue = sample_name.split('_rep')[0]
        parquet_path = os.path.join(input_dir_base, tissue, f"{sample_name}.parquet")
        
        print(f"  - Processing {sample_name}...")
        
        if not os.path.exists(parquet_path):
            print(f"    - Warning: Input file not found at {parquet_path}. Skipping.")
            continue
            
        sample_df = pd.read_parquet(parquet_path)
        
        plot_stacked_contribution(sample_df, name=sample_name, output_dir=output_dir)
        plot_bar_genecount(sample_df, name=sample_name, output_dir=output_dir)
        
    print(f"\nAnalysis complete. Plots saved to {output_dir}")

if __name__ == "__main__":
    main()