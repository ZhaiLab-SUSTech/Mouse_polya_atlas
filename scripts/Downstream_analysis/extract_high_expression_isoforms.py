# -*- coding: utf-8 -*-
import os
import sys
from typing import List

import pandas as pd

def get_high_expression_isoforms(
    info_df_path: str, 
    min_sum: int, 
    min_tissues: int
) -> List[str]:
    """
    Filters the main info DataFrame to identify and return a list of
    highly and widely expressed isoforms.
    """
    print(f"Loading merged info DataFrame from {info_df_path}...")
    info_df = pd.read_pickle(info_df_path)
    
    print(f"Filtering for isoforms with expression sum > {min_sum} in > {min_tissues} tissues...")
    gene_stats = info_df.groupby('gene_info').agg(
        expression_sum=('expression', 'sum'),
        sample_num=('tissue', 'nunique')
    )
    
    highly_expressed_mask = (
        (gene_stats['expression_sum'] > min_sum) &
        (gene_stats['sample_num'] > min_tissues) &
        (gene_stats.index.str.contains("\.", na=False))
    )
    
    high_expression_isoforms = gene_stats[highly_expressed_mask].index.tolist()
    print(f"Found {len(high_expression_isoforms)} isoforms passing the filters.")
    return high_expression_isoforms

def main():
    """
    Main function to generate and save the list of high-expression isoforms.
    """
    # --- Configuration ---
    merged_info_path = "results/06_merged_data/merged_info_df.pkl"
    output_dir = "results/downstream_analysis/"
    output_file_path = os.path.join(output_dir, "high_expression_isoforms.txt")

    # Parameters for filtering, consistent with previous scripts
    min_expression_sum = 240
    min_tissue_count = 16

    # --- Workflow ---
    os.makedirs(output_dir, exist_ok=True)
    
    isoform_list = get_high_expression_isoforms(merged_info_path, min_expression_sum, min_tissue_count)
    
    print(f"Saving list to {output_file_path}...")
    with open(output_file_path, 'w') as f:
        for isoform_id in isoform_list:
            f.write(f"{isoform_id}\n")
            

if __name__ == "__main__":
    main()