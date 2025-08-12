#!/bin-bash
#PBS -N calc_distance_matrices
#PBS -l nodes=1:ppn=36,mem=120G
#PBS -j oe
#PBS -l walltime=200:00:00

set -euo pipefail
cd "$PBS_O_WORKDIR"

# --- Configuration ---
THREADS=36

# Input files
MERGED_INFO_DF="results/06_merged_data/merged_info_df.pkl"
MERGED_ARRAY="results/06_merged_data/merged_array.npy"

COMBINED_CPM="results/03_polya_matrix/combined_cpm_polya_matrix.csv"

# Scripts 
POLYA_COPH_SCRIPT="scripts/preprocessing_pipeline/07_distance_matrix/calculate_polya_cophenet_distance.py"
POLYA_PEARSON_SCRIPT="scripts/preprocessing_pipeline/07_distance_matrix/calculate_polya_pearson_distance.py"
MERGE_SCRIPT="scripts/preprocessing_pipeline/07_distance_matrix/merge_distance.py"

EXPR_SCRIPT="scripts/preprocessing_pipeline/07_distance_matrix/calculate_expression_distance.py"

# Output files
POLYA_COPH_H5="results/07_distance_matrix/polya_cophenetic_pearson.h5"
POLYA_PEARSON_H5="results/07_distance_matrix/polya_direct_pearson.h5"
POLYA_FINAL_H5="results/07_distance_matrix/polya_distance.h5"

EXPR_FINAL_H5="results/07_distance_matrix/expression_distance.h5"

# --- Execution ---
python3 "$POLYA_COPH_SCRIPT" --info_df "$MERGED_INFO_DF" --array "$MERGED_ARRAY" --output_h5 "$POLYA_COPH_H5" --threads "$THREADS"

python3 "$POLYA_PEARSON_SCRIPT" --info_df "$MERGED_INFO_DF" --array "$MERGED_ARRAY" --output_h5 "$POLYA_PEARSON_H5" --threads "$THREADS"

python3 "$MERGE_SCRIPT" --matrix1_h5 "$POLYA_COPH_H5" --matrix2_h5 "$POLYA_PEARSON_H5" --output_h5 "$POLYA_FINAL_H5"

python3 "$EXPR_SCRIPT" --info_df "$MERGED_INFO_DF" --cpm_matrix "$COMBINED_CPM" --output_h5 "$EXPR_FINAL_H5" --threads "$THREADS"

echo "All distance matrix calculations for Step 07 are complete."