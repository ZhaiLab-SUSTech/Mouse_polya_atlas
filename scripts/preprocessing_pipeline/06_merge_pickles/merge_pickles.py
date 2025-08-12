# -*- coding: utf-8 -*-
import os
import sys
import pickle
import argparse
from typing import List, Dict

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'Downstream_analysis')))
from config_downstream import TISSUE_LIST

def load_tissue_data(pkl_dir: str, tissue: str) -> Dict:
    base_path = os.path.join(pkl_dir, tissue, f"{tissue}")
    try:
        with open(f'{base_path}_dataarray.pkl', 'rb') as f:
            data = pickle.load(f)
        with open(f'{base_path}_sampleinfo.pkl', 'rb') as f:
            sample_info = pickle.load(f)
        with open(f'{base_path}_geneinfo_isoform.pkl', 'rb') as f:
            gene_info = pickle.load(f)
        with open(f'{base_path}_expression.pkl', 'rb') as f:
            expression = pickle.load(f)
        with open(f'{base_path}_geneinfo_geneid.pkl', 'rb') as f:
            gene_id = pickle.load(f)
        with open(f'{base_path}_geneinfo_genesymbol.pkl', 'rb') as f:
            gene_symbol = pickle.load(f)
        
        return {
            "data": data, "sample_info": sample_info, "gene_info": gene_info,
            "expression": expression, "gene_id": gene_id, "gene_symbol": gene_symbol
        }
    except FileNotFoundError as e:
        print(f"Warning: Could not load all files for tissue '{tissue}'. Skipping. Details: {e}", file=sys.stderr)
        return None

def create_info_dataframe(tissue_data: Dict, tissue_name: str) -> pd.DataFrame:
    df = pd.DataFrame({
        "sample_info": tissue_data["sample_info"],
        "gene_info": tissue_data["gene_info"],
        "expression": tissue_data["expression"],
        "gene_id": tissue_data["gene_id"],
        "gene_symbol": tissue_data["gene_symbol"],
        "tissue": tissue_name
    })
    return df

def main(input_dir: str, output_info_path: str, output_array_path: str) -> None:
    all_info_dfs: List[pd.DataFrame] = []
    all_arrays: List[np.ndarray] = []

    for tissue in TISSUE_LIST:
        print(f"  - Processing {tissue}...")
        tissue_data = load_tissue_data(input_dir, tissue)
        
        if tissue_data:
            # Create info DataFrame and append to list
            info_df = create_info_dataframe(tissue_data, tissue)
            all_info_dfs.append(info_df)
            
            # The data is a list of arrays, concatenate them first
            concatenated_array = np.concatenate(tissue_data["data"], axis=0)
            all_arrays.append(concatenated_array)

    # Final merge after the loop
    merged_info_df = pd.concat(all_info_dfs, ignore_index=True)
    merged_array = np.concatenate(all_arrays, axis=0)
    
    # Filter out rows with invalid gene information
    valid_rows_mask = merged_info_df['gene_info'].str.startswith("EN", na=False)
    final_info_df = merged_info_df[valid_rows_mask].reset_index(drop=True)
    final_array = merged_array[valid_rows_mask]

    # Create output directories and save files
    os.makedirs(os.path.dirname(output_info_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_array_path), exist_ok=True)
    
    print(f"Saving merged info DataFrame to {output_info_path}...")
    final_info_df.to_pickle(output_info_path)
    
    print(f"Saving merged numpy array to {output_array_path}...")
    np.save(output_array_path, final_array)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge data from tissue-specific pickle files.")
    parser.add_argument('--input_dir', required=True, help="Directory containing the output of the 'extract_polya' step.")
    parser.add_argument('--output_info', required=True, help="Path to save the merged info DataFrame (as .pkl).")
    parser.add_argument('--output_array', required=True, help="Path to save the merged numpy array (as .npy).")
    args = parser.parse_args()
    main(args.input_dir, args.output_info, args.output_array)