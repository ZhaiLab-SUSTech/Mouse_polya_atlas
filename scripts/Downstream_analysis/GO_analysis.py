# -*- coding: utf-8 -*-
import os
from typing import List, Dict, Optional

import numpy as np
import pandas as pd
import gseapy as gp
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns

def set_plotting_style() -> None:
    sns.set(font='Arial')
    plt.rcParams['svg.fonttype'] = 'none'
    style = sns.axes_style('white')
    style.update(sns.axes_style('ticks'))
    sns.set(font_scale=1.4, style=style)

def run_gsea_enrichment(
    modules_df: pd.DataFrame, 
    gmt_path: str,
    module_col_name: str
) -> pd.DataFrame:
    enrichment_results: Dict[str, pd.DataFrame] = {}
    
    # Dynamically find all unique modules to analyze
    unique_modules = sorted(modules_df[module_col_name].unique())
    
    for module_id in unique_modules:
        
        gene_list = modules_df[modules_df[module_col_name] == module_id]['gene_symbol'].tolist()
        
        if not gene_list:
            print(f"  - Module {module_id}: Skipping, no genes found.")
            continue
            
        try:
            print(f"  - Analyzing Module {module_id} ({len(gene_list)} genes)...")
            enr_res = gp.enrichr(
                gene_list=gene_list,
                gene_sets=[gmt_path],
                organism='Mouse',
            )
            if not enr_res.results.empty:
                enrichment_results[f'module_{module_id}'] = enr_res.results
        except Exception as e:
            print(f"  - Module {module_id}: Could not be processed. Error: {e}")

    all_modules_list = []
    for module_name, df in enrichment_results.items():
        df[['n_intersect', 'n_total']] = df['Overlap'].str.split('/', expand=True).astype(int)
        df['Gene Ratio'] = df['n_intersect'] / df['n_total']
        df['-log10(padj)'] = -np.log10(df['Adjusted P-value'])
        
        # Filter for significant terms and take the top 10
        df_filtered = df[df['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value').head(10)
        
        df_filtered['Module'] = module_name.split('_')[1]
        all_modules_list.append(df_filtered)
        
    if not all_modules_list:
        print("Warning: No significant enrichment results found across all modules.")
        return pd.DataFrame()

    return pd.concat(all_modules_list, ignore_index=True)

def cluster_and_order_terms(combined_df: pd.DataFrame) -> pd.DataFrame:
    # Filter for terms to be plotted
    plot_df = combined_df[
        (combined_df['Adjusted P-value'] < 0.05) & 
        (combined_df['n_intersect'] >= 3)
    ].copy()
    
    if plot_df.empty:
        print("Warning: No terms left after filtering for clustering. Cannot create plot.")
        return pd.DataFrame()

    term_module_matrix = pd.crosstab(plot_df['Term'], plot_df['Module']).astype(bool).astype(int)
    
    distance_matrix = pdist(term_module_matrix.values, metric='jaccard')
    Z = linkage(distance_matrix, method='average')
    term_order = term_module_matrix.index[leaves_list(Z)]
    
    plot_df['Term'] = pd.Categorical(plot_df['Term'], categories=term_order, ordered=True)
    plot_df.sort_values(['Term', 'Module'], inplace=True)
    
    return plot_df

def plot_go_dotplot(
    plot_df: pd.DataFrame, 
    output_path: str,
    title: str,
    module_order: Optional[List[str]] = None
) -> None:
    # Scale point sizes
    s_min, s_max = 20, 200
    p_min, p_max = plot_df['Gene Ratio'].min(), plot_df['Gene Ratio'].max()
    if p_max - p_min == 0:
        plot_df['scaled_size'] = s_max
    else:
        plot_df['scaled_size'] = s_min + (plot_df['Gene Ratio'] - p_min) / (p_max - p_min) * (s_max - s_min)

    plt.figure(figsize=(12, max(8, len(plot_df['Term'].unique()) * 0.25)))
    ax = plt.gca()

    if module_order is None:
        module_order = sorted(plot_df['Module'].unique(), key=int)
    
    x_positions = {str(m): i for i, m in enumerate(module_order)}
    plot_df['x_pos'] = plot_df['Module'].map(x_positions)

    scatter = ax.scatter(
        x=plot_df['x_pos'], y=plot_df['Term'], s=plot_df['scaled_size'],
        c=plot_df['-log10(padj)'], cmap='RdYlBu_r', alpha=0.8,
        edgecolor='black', linewidth=0.5, vmin=1, vmax=10
    )
    
    ax.set_xticks(range(len(module_order)))
    ax.set_xticklabels([f'Module {m}' for m in module_order], rotation=45, ha='right')
    ax.set_xlabel('Module', fontsize=12)
    ax.set_ylabel('GO Term', fontsize=12)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    cbar = plt.colorbar(scatter, pad=0.02, aspect=30)
    cbar.set_label('-log10(Adjusted P-value)', fontsize=12)

    legend_sizes_raw = np.linspace(p_min, p_max, 3)
    legend_sizes_scaled = np.linspace(s_min, s_max, 3)

    legend_handles = []
    for i, raw_size in enumerate(legend_sizes_raw):
        handle = ax.scatter([], [], s=legend_sizes_scaled[i],
                            label=f'{raw_size:.2f}', # Format label to 2 decimal places
                            color='gray', alpha=0.7,
                            edgecolor='black', linewidth=0.5)
        legend_handles.append(handle)

    size_legend = ax.legend(
        handles=legend_handles,
        title='Gene Ratio',
        loc='upper left',
        bbox_to_anchor=(1.1, 0.5), # Position legend outside the plot
        labelspacing=1.5,
        borderpad=1,
        handletextpad=1.5,
        fontsize=10
    )
    size_legend.get_title().set_fontsize('12')

    plt.title(title, fontsize=16)
    plt.subplots_adjust(right=0.85)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {output_path}")