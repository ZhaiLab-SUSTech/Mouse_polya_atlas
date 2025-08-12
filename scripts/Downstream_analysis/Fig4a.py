# -*- coding: utf-8 -*-
import os
import sys
import subprocess

import numpy as np
import pandas as pd
import h5py
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import marsilea as ma
import marsilea.plotter as mp

def run_r_script(r_script_path: str, input_csv: str, output_dir: str) -> bool:
    rscript_exe = "Rscript"
    if not any(os.access(os.path.join(path, rscript_exe), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
        print(f"Error: '{rscript_exe}' not found in PATH.", file=sys.stderr)
        return False

    command = [rscript_exe, r_script_path, input_csv, output_dir]
    print(f"Executing command: {' '.join(command)}")
    
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print("R script executed successfully.")
        return True
    except FileNotFoundError:
        print(f"Error: R script not found at {r_script_path}", file=sys.stderr)
        return False
    except subprocess.CalledProcessError as e:
        print(f"Error during R script execution: {e}", file=sys.stderr)
        print(f"R stderr:\n{e.stderr}")
        return False

def plot_final_heatmap(
    submatrix_path: str, 
    modules_path: str, 
    output_path: str
) -> None:
    reordered_matrix_df = pd.read_csv(submatrix_path, index_col=0)
    wgcna_module = pd.read_csv(modules_path)

    gene_cat = wgcna_module['Module'].tolist()
    cat_order = sorted(wgcna_module['Module'].unique())
    
    colors = sns.color_palette("husl", n_colors=len(cat_order))
    hex_colors = [mcolors.rgb2hex(color) for color in colors]
    
    h = ma.Heatmap(
        reordered_matrix_df.values, cmap="viridis", width=6, height=6, vmin=0, vmax=1.8
    )
    
    h.hsplit(labels=gene_cat, order=cat_order, spacing=0.005)
    h.vsplit(labels=gene_cat, order=cat_order, spacing=0.005)
    
    h.add_left(mp.Chunk(cat_order, colors=hex_colors, fontsize=5), pad=0.05)
    h.add_top(mp.Chunk([""]*len(cat_order), colors=hex_colors, fontsize=5), pad=0.05)
    
    h.add_dendrogram("left", colors=hex_colors, method='single', metric='correlation')
    h.add_dendrogram("top", colors=hex_colors, method='single', metric='correlation')

    h.render()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Final heatmap saved to {output_path}")

def main():
    merged_info_path = "results/06_merged_data/merged_info_df.pkl"
    distance_matrix_path = "results/07_distance_matrix/pearson_expression_acrosstissue.h5"
    
    python_output_dir = "results/downstream_analysis/07_submatrix"
    submatrix_csv_path = os.path.join(python_output_dir, "clustered_expression_submatrix.csv")
    
    r_script_path = "scripts/Downstream_analysis/cutree_expression.R"
    r_output_dir = "results/downstream_analysis/06_expression_modules"
    modules_csv_path = os.path.join(r_output_dir, "expression_modules.csv")

    final_heatmap_path = os.path.join(r_output_dir, "final_expression_annotated_heatmap.pdf")
    high_expression_isoforms_path = "results/downstream_analysis/high_expression_isoforms.txt"

    info_df = pd.read_pickle(merged_info_path)
    with h5py.File(distance_matrix_path, 'r') as f:
        expression_matrix = f['distance_matrix'][:]

    isoform_map = info_df.drop_duplicates('gene_info').set_index('gene_info')
    isoform2symbol = isoform_map['gene_symbol'].to_dict()

    if not os.path.exists(high_expression_isoforms_path):
        print(f"Error: High expression isoform list not found.", file=sys.stderr)
        print(f"Please run 'extract_high_expression_isoforms.py' first.", file=sys.stderr)
        sys.exit(1)

    with open(high_expression_isoforms_path, 'r') as f:
        high_expression_genes = [line.strip() for line in f]
    
    all_unique_genes = info_df['gene_info'].unique().tolist()
    gene_to_order_map = {gene: i for i, gene in enumerate(all_unique_genes)}
    high_expression_idx = [gene_to_order_map[g] for g in high_expression_genes if g in gene_to_order_map]
    
    if len(high_expression_idx) < 2: return
        
    sub_matrix = expression_matrix[np.ix_(high_expression_idx, high_expression_idx)]
    Z = linkage(sub_matrix, method='average', metric='euclidean')
    leaves_order = leaves_list(Z)
    
    reordered_matrix = sub_matrix[leaves_order, :][:, leaves_order]
    reordered_isoforms = [high_expression_genes[i] for i in leaves_order]
    reordered_symbols = [isoform2symbol.get(iso, '') for iso in reordered_isoforms]
    final_index_names = [f"{sym}-{iso}" for sym, iso in zip(reordered_symbols, reordered_isoforms)]
    
    expression_sub_df = pd.DataFrame(reordered_matrix, index=final_index_names, columns=final_index_names)
    os.makedirs(python_output_dir, exist_ok=True)
    expression_sub_df.to_csv(submatrix_csv_path)
    print(f"Clustered expression sub-matrix saved to {submatrix_csv_path}")

    r_success = run_r_script(r_script_path, submatrix_csv_path, r_output_dir)

    if r_success and os.path.exists(modules_csv_path):
        plot_final_heatmap(
            submatrix_path=submatrix_csv_path,
            modules_path=modules_csv_path,
            output_path=final_heatmap_path
        )
    else:
        print("Skipping final heatmap generation due to errors or missing files.", file=sys.stderr)


if __name__ == "__main__":
    main()