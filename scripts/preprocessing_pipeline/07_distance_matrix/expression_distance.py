# -*- coding: utf-8 -*-
import argparse
from typing import Tuple

import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

def generate_cpm_vector(isoform_id: str, cpm_matrix: pd.DataFrame) -> np.ndarray:
    short_id = isoform_id.split('.')[0]
    try:
        return cpm_matrix.loc[short_id].values
    except KeyError:
        return np.zeros(cpm_matrix.shape[1])

def process_full_row(i: int, matrix: np.ndarray) -> Tuple[int, np.ndarray]:
    corr = np.corrcoef(matrix[i], matrix)[0, :]
    dist = np.sqrt(2 * (1 - np.clip(corr, -1.0, 1.0))).astype(np.float32)
    dist[i] = 0.0
    return i, dist

def main(args):
    info_df = pd.read_pickle(args.info_df)
    cpm_matrix = pd.read_csv(args.cpm_matrix)

    # Create a short isoform ID (without version) for indexing
    cpm_matrix['new_isoform_id'] = cpm_matrix['isoform_id'].str.split('.').str[0]
    cpm_matrix.set_index('new_isoform_id', inplace=True)
    cpm_matrix_trim = cpm_matrix.filter(regex='_bam_cpm$')
    
    isoform_list = info_df.gene_info.unique().tolist()

    print(f"Extracting CPM vectors for {len(isoform_list)} isoforms...")
    cpm_vectors_matrix = np.array(Parallel(n_jobs=args.threads)(
        delayed(generate_cpm_vector)(iso, cpm_matrix_trim) for iso in isoform_list
    ), dtype=np.float32)

    print("Computing final distance matrix between all isoforms...")
    num_isoforms = cpm_vectors_matrix.shape[0]
    with h5py.File(args.output_h5, 'w') as hf:
        dset = hf.create_dataset('distance_matrix', (num_isoforms, num_isoforms), 'f4', chunks=(args.block_size, args.block_size), compression='gzip')
        for start in range(0, num_isoforms, args.block_size):
            end = min(start + args.block_size, num_isoforms)
            results = Parallel(n_jobs=args.threads)(delayed(process_full_row)(i, cpm_vectors_matrix) for i in range(start, end))
            for i, dist_row in results:
                dset[i, :] = dist_row
            hf.flush()
            
    print(f"Expression distance matrix saved to {args.output_h5}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate expression-based distance matrix from CPM values.")
    parser.add_argument('--info_df', required=True, help="Path to merged info DataFrame.")
    parser.add_argument('--cpm_matrix', required=True, help="Path to combined CPM matrix CSV.")
    parser.add_argument('--output_h5', required=True, help="Path to save the output H5 distance matrix.")
    parser.add_argument('--block_size', type=int, default=2000)
    parser.add_argument('--threads', type=int, default=8)
    main(parser.parse_args())