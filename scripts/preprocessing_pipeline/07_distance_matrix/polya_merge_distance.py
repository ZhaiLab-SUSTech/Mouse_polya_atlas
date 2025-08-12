# -*- coding: utf-8 -*-
import argparse
import h5py
import numpy as np

def main(args):
    """
    Merges two HDF5 distance matrices block by block using the fmin function.
    """
    print(f"Merging matrices:\n  - {args.matrix1_h5}\n  - {args.matrix2_h5}")
    
    with h5py.File(args.matrix1_h5, 'r') as hf1, \
         h5py.File(args.matrix2_h5, 'r') as hf2, \
         h5py.File(args.output_h5, 'w') as hf_merged:

        dset1 = hf1['distance_matrix']
        dset2 = hf2['distance_matrix']
        
        if dset1.shape != dset2.shape:
            raise ValueError("Input distance matrices must have the same shape.")
            
        num_genes = dset1.shape[0]

        dset_merged = hf_merged.create_dataset(
            'distance_matrix', 
            shape=(num_genes, num_genes), 
            dtype='float32', 
            chunks=(args.block_size, args.block_size), 
            compression='gzip'
        )

        for start in range(0, num_genes, args.block_size):
            end = min(start + args.block_size, num_genes)
            print(f"  - Merging block {start}-{end}...")
            
            block1 = dset1[start:end, :]
            block2 = dset2[start:end, :]
            
            # Use np.fmin to merge the blocks, as per the original script's logic
            merged_block = np.fmin(block1, block2)
            dset_merged[start:end, :] = merged_block
            hf_merged.flush()
            
    print(f"Merged distance matrix has been saved to {args.output_h5}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge two HDF5 distance matrices.")
    parser.add_argument('--matrix1_h5', required=True, help="Path to the first HDF5 distance matrix.")
    parser.add_argument('--matrix2_h5', required=True, help="Path to the second HDF5 distance matrix.")
    parser.add_argument('--output_h5', required=True, help="Path to save the merged HDF5 distance matrix.")
    parser.add_argument('--block_size', type=int, default=2000)
    main(parser.parse_args())