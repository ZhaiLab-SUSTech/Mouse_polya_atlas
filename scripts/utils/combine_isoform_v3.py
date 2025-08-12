#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import logging
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from merge_utils import optimize_dataframe, parse_gtf, read_read_info

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(filename)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)


def process_sample(args):
    rep_name, results_dir, isoquant_dir, gtf_dir, output_dir = args
    try:
        tissue, rep = rep_name.rsplit("_", 1)
    except ValueError:
        print(
            f"Error: Sample name {rep_name} does not follow the expected format 'tissue_rep'.",
            file=sys.stderr,
        )
        return None
    print(f"Processing sample: {rep_name}", file=sys.stderr)

    sample_iso_path = os.path.join(
        isoquant_dir, rep_name, rep_name, f"{rep_name}.read_assignments.tsv.gz"
    )
    if not os.path.exists(sample_iso_path):
        print(f"Error: {sample_iso_path} does not exist.", file=sys.stderr)
        return None

    try:
        with gzip.open(sample_iso_path, "rt") as f:
            sample_iso = pd.read_csv(
                f, sep="\t", skiprows=2, dtype=str
            )  
            expected_columns = [
                "read_id",
                "chr",
                "strand",
                "isoform_id",
                "gene_id",
                "assignment_type",
                "assignment_events",
                "exons",
                "additional_info",
            ]
            if len(sample_iso.columns) < len(expected_columns):
                print(
                    f"Error: Unexpected number of columns in {sample_iso_path}. Expected at least {len(expected_columns)}, got {len(sample_iso.columns)}.",
                    file=sys.stderr,
                )
                return None
            sample_iso.columns = expected_columns[: len(sample_iso.columns)]
            logging.info("load isoquant file")
    except Exception as e:
        print(f"Error reading {sample_iso_path}: {e}", file=sys.stderr)
        return None

    if "read_id" not in sample_iso.columns or "chr" not in sample_iso.columns:
        print(
            f"Error: 'read_id' or 'chr' column missing in {sample_iso_path}.",
            file=sys.stderr,
        )
        return None

    gtf_path = os.path.join(gtf_dir, f"{rep_name}.extended_annotation.gtf")
    if not os.path.exists(gtf_path):
        print(f"Error: {gtf_path} does not exist.", file=sys.stderr)
        return None

    try:
        _, transcript_to_gene, gene_to_symbol = parse_gtf(gtf_path)
        logging.info("load gtf")
    except Exception as e:
        print(f"Error parsing GTF file {gtf_path}: {e}", file=sys.stderr)
        return None

    read_info_path = os.path.join(results_dir, f"{rep_name}.read.info.txt")
    if not os.path.exists(read_info_path):
        print(f"Error: {read_info_path} does not exist.", file=sys.stderr)
        return None

    read_info_df = read_read_info(read_info_path)
    logging.info("load read info")

    read_info_df["read_id"] = read_info_df["read_core_id"].str.split(",").str[0]

    read_info_df["read_id"] = read_info_df["read_id"].astype(str)
    read_info_df["chr"] = read_info_df["chr"].astype(str)
    sample_iso["read_id"] = sample_iso["read_id"].astype(str)
    sample_iso["chr"] = sample_iso["chr"].astype(str)

    merged = pd.merge(sample_iso, read_info_df, on=["read_id", "chr"], how="outer")
    logging.info("merged")

    if merged.empty:
        print(f"Warning: No data after merging for sample {rep_name}.", file=sys.stderr)
        return None

    try:
        merged["gene_symbol"] = merged["gene_id"].apply(
            lambda x: gene_to_symbol.get(x, np.nan)
        )
    except Exception as e:
        print(f"Error adding gene_symbol for sample {rep_name}: {e}", file=sys.stderr)
        merged["gene_symbol"] = np.nan

    try:
        merged["exons_split"] = merged["exons"].str.split("-", n=1000, expand=False)
    except Exception as e:
        print(f"Error splitting exons for sample {rep_name}: {e}", file=sys.stderr)
        merged["exons_split"] = [[] for _ in range(len(merged))]

    try:
        merged["exons_ordered"] = merged.apply(
            lambda row: (
                row["exons_split"][::-1] if row["strand"] == "-" else row["exons_split"]
            ),
            axis=1,
        )
    except Exception as e:
        print(f"Error ordering exons for sample {rep_name}: {e}", file=sys.stderr)
        merged["exons_ordered"] = merged["exons_split"]

    try:
        merged["end_site"] = merged["exons_ordered"].apply(
            lambda x: x[-1] if isinstance(x, list) and len(x) > 0 else np.nan
        )
    except Exception as e:
        print(f"Error extracting end_site for sample {rep_name}: {e}", file=sys.stderr)
        merged["end_site"] = np.nan

    try:
        merged["end_site"] = pd.to_numeric(merged["end_site"], errors="coerce")
    except Exception as e:
        print(
            f"Error converting end_site to numeric for sample {rep_name}: {e}",
            file=sys.stderr,
        )
        merged["end_site"] = np.nan

    merged = merged.drop(["exons_ordered", "exons_split"], axis=1, errors="ignore")

    logging.info("Optimizing df...")

    try:
        merged = optimize_dataframe(merged)
    except Exception as e:
        print(f"Error optimizing DataFrame for sample {rep_name}: {e}", file=sys.stderr)

    for col in merged.columns:
        if merged[col].dtype == "object":
            merged[col] = merged[col].astype(str)
        elif pd.api.types.is_integer_dtype(merged[col]):
            merged[col] = merged[col].astype("Int64") 
        elif pd.api.types.is_float_dtype(merged[col]):
            merged[col] = merged[col].astype("float64")
        elif pd.api.types.is_categorical_dtype(merged[col]):
            pass 
        else:
            pass

    logging.info("finish process, start saving")

    if output_dir:
        try:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"{rep_name}.parquet")
            merged.to_parquet(output_path, index=False, compression="snappy")
            print(f"Saved processed data to {output_path}", file=sys.stderr)
        except Exception as e:
            print(f"Error saving Parquet for sample {rep_name}: {e}", file=sys.stderr)
            return None



def main():
    parser = argparse.ArgumentParser(description="Merge and process samples.")
    parser.add_argument(
        "--replicates",
        nargs="+",
        required=True,
        help="List of replicate sample names, e.g., tissue_rep1 tissue_rep2",
    )
    parser.add_argument(
        "--results_dir", required=True, help="Path to results directory"
    )
    parser.add_argument(
        "--isoquant_dir", required=True, help="Path to Isoquant_raw directory"
    )
    parser.add_argument(
        "--gtf_dir", required=True, help="Path to discovery_gtf directory"
    )
    parser.add_argument("--output_dir", required=True, help="Path to save merged data")
    parser.add_argument(
        "--num_processes", type=int, default=1, help="Number of parallel processes"
    )

    args = parser.parse_args()

    replicates = args.replicates
    results_dir = args.results_dir
    isoquant_dir = args.isoquant_dir
    gtf_dir = args.gtf_dir
    output_dir = args.output_dir
    num_processes = args.num_processes

    pool_args = [
        (rep_name, results_dir, isoquant_dir, gtf_dir, output_dir)
        for rep_name in replicates
    ]

    with Pool(processes=num_processes) as pool:
        pool.map(process_sample, pool_args)

if __name__ == "__main__":
    main()
