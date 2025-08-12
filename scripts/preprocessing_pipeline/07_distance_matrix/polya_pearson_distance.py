# -*- coding: utf-8 -*-
import argparse
from typing import Tuple

import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from distance_utils import generate_gene_array

def generate_gene_array_flat(gene: str, info_df: pd.DataFrame, data_array: np.ndarray) -> np.ndarray:
    result_array = generate_gene_array(gene, info_df, data_array)
    return np.nan_to_num(result_array.flatten(), nan=0)

def process_full_row(i: int, matrix: np.ndarray) -> Tuple[int, np.ndarray]:
    corr = np.corrcoef(matrix[i], matrix)[0, :]
    dist = np.sqrt(2 * (1 - np.clip(corr, -1.0, 1.0))).astype(np.float32)
    dist[i] = 0.0
    return i, dist

def main(args):
    info_df = pd.read_pickle(args.info_df)
    data_array = np.load(args.array)
    gene_list = info_df.gene_info.unique().tolist()

    flattened_matrix = np.array(Parallel(n_jobs=args.threads)(
        delayed(generate_gene_array_flat)(g, info_df, data_array) for g in gene_list
    ), dtype=np.float32)

    num_genes = flattened_matrix.shape[0]
    with h5py.File(args.output_h5, 'w') as hf:
        dset = hf.create_dataset('distance_matrix', (num_genes, num_genes), 'f4', chunks=(args.block_size, args.block_size), compression='gzip')
        for start in range(0, num_genes, args.block_size):
            end = min(start + args.block_size, num_genes)
            print(f"  - 处理块 {start}-{end}...")
            results = Parallel(n_jobs=args.threads)(delayed(process_full_row)(i, flattened_matrix) for i in range(start, end))
            for i, dist_row in results:
                dset[i, :] = dist_row
            hf.flush()
            
    print(f"PolyA pearson distance matrix saved to {args.output_h5}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate direct Pearson distance on Poly(A) distributions.")
    parser.add_argument('--info_df', required=True)
    parser.add_argument('--array', required=True)
    parser.add_argument('--output_h5', required=True)
    parser.add_argument('--block_size', type=int, default=2000)
    parser.add_argument('--threads', type=int, default=8)
    main(parser.parse_args())