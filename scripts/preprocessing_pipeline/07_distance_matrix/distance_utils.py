# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'Downstream_analysis')))
from config_downstream import TISSUE_LIST

def generate_gene_array(gene: str, info_df: pd.DataFrame, data_array: np.ndarray) -> np.ndarray:
    gene_df = info_df[info_df['gene_info'] == gene]
    # Create an empty array with NaNs to fill
    result_array = np.full((len(TISSUE_LIST), data_array.shape[1]), np.nan)
    
    # Create a mapping of tissue to its index in the original data_array
    tissue_to_pos = {tissue: idx for idx, tissue in zip(gene_df.index, gene_df.sample_info)}
    
    # Fill the result_array based on the tissue mapping
    for i, tissue in enumerate(TISSUE_LIST):
        pos = tissue_to_pos.get(tissue)
        if pos is not None:
            result_array[i] = data_array[pos]
            
    return result_array