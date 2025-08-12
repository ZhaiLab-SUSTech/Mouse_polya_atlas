# -*- coding: utf-8 -*-
import argparse
from itertools import combinations
from typing import Tuple

import h5py
import numpy as np
import pandas as pd
import ot
from joblib import Parallel, delayed
from scipy.cluster.hierarchy import linkage, cophenet

from distance_utils import generate_gene_array

def compute_wasserstein_dist(array: np.ndarray, bin_width: int) -> np.ndarray:
    num_tissues = array.shape[0]
    dist_vector = []
    bins = np.arange(10, 301, bin_width, dtype=np.float64)
    for i, j in combinations(range(num_tissues), 2):
        u, v = array[i], array[j]
        u_sum, v_sum = u.sum(), v.sum()
        if u_sum == 0 and v_sum == 0: dist = 0.0
        elif u_sum == 0 or v_sum == 0: dist = np.inf
        else: dist = ot.wasserstein_1d(bins, bins, u / u_sum, v / v_sum, p=2)
        dist_vector.append(dist)
    dist_vector = np.array(dist_vector, dtype=np.float32)
    if np.any(np.isinf(dist_vector)):
        max_finite = np.max(dist_vector[np.isfinite(dist_vector)], initial=1e6)
        dist_vector[np.isinf(dist_vector)] = max_finite * 1.2
    return dist_vector

def compute_cophenet_vector(array: np.ndarray, bin_width: int) -> np.ndarray:
    array = np.nan_to_num(array, nan=0)
    condensed_dist = compute_wasserstein_dist(array, bin_width)
    Z = linkage(condensed_dist, method='complete')
    _, coph_dist = cophenet(Z, condensed_dist)
    return coph_dist.astype(np.float32)

def process_full_row(i: int, matrix: np.ndarray) -> Tuple[int, np.ndarray]:
    corr = np.corrcoef(matrix[i], matrix)[0, :]
    dist = np.sqrt(2 * (1 - np.clip(corr, -1.0, 1.0))).astype(np.float32)
    dist[i] = 0.0
    return i, dist

def main(args):
    info_df = pd.read_pickle(args.info_df)
    data_array = np.load(args.array)
    gene_list = info_df.gene_info.unique().tolist()
    
    print("Computing cophenetic vectors (PolyA)...")
    cophenet_matrix = np.array(Parallel(n_jobs=args.threads)(
        delayed(compute_cophenet_vector)(generate_gene_array(g, info_df, data_array), args.bin_width)
        for g in gene_list
    ), dtype=np.float32)

    print("Computing final distance matrix (PolyA)...")
    num_genes = cophenet_matrix.shape[0]
    with h5py.File(args.output_h5, 'w') as hf:
        dset = hf.create_dataset('distance_matrix', (num_genes, num_genes), 'f4', chunks=(args.block_size, args.block_size), compression='gzip')
        for start in range(0, num_genes, args.block_size):
            end = min(start + args.block_size, num_genes)
            results = Parallel(n_jobs=args.threads)(delayed(process_full_row)(i, cophenet_matrix) for i in range(start, end))
            for i, dist_row in results:
                dset[i, :] = dist_row
            hf.flush()
    print(f"PolyA distance matrix saved to {args.output_h5}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate Poly(A)-based distance matrix.")
    parser.add_argument('--info_df', required=True)
    parser.add_argument('--array', required=True)
    parser.add_argument('--output_h5', required=True)
    parser.add_argument('--bin_width', type=int, default=5)
    parser.add_argument('--block_size', type=int, default=2000)
    parser.add_argument('--threads', type=int, default=8)
    main(parser.parse_args())