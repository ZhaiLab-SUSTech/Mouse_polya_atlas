# -*- coding: utf-8 -*-
import argparse
import functools
import pandas as pd
from glob import glob
import os

def read_and_order_cols(filename: str) -> pd.DataFrame:
    """Reads a CPM file and ensures required columns are present."""
    df = pd.read_csv(filename, sep="\t")
    required_cols = ["isoform_id", "gene_id", "gene_symbol"]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"File {filename} is missing required columns.")
    return df

def combine_cpm_files(input_dir: str, output_file: str):
    """Finds all CPM files in a directory, merges them, and saves the result."""
    # Find all per-sample/per-tissue cpm files
    all_files = glob(os.path.join(input_dir, "*/*_cpm_polya_gene.csv"))
    if not all_files:
        raise FileNotFoundError(f"No '*_cpm_polya_gene.csv' files found in subdirectories of {input_dir}")
    
    print(f"Found {len(all_files)} CPM files to merge.")
    
    dfs = [read_and_order_cols(f) for f in all_files]
    
    # Merge all dataframes on the key columns
    combined_df = functools.reduce(
        lambda left, right: pd.merge(
            left, right, on=["isoform_id", "gene_id", "gene_symbol"], how="outer"
        ),
        dfs,
    ).fillna(0)
    
    print(f"Saving combined CPM matrix to {output_file}")
    combined_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine multiple CPM matrix files.")
    parser.add_argument('--input_dir', required=True, help="Directory containing the output of step 03 (e.g., 'results/03_polya_matrix').")
    parser.add_argument('--output_file', required=True, help="Path to save the combined CSV file.")
    args = parser.parse_args()
    combine_cpm_files(args.input_dir, args.output_file)