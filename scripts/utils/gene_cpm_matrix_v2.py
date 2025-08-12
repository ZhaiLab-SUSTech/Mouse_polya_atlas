"""
Description: calcuate polya_matrix; filter full_length and pa>10; 
"""

import argparse
import functools
import logging
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from merge_utils import parse_gtf, unify_columns
from sample_parse_utils import parse_sample

def configure_logging(log_file_path: str):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file_path),
            logging.StreamHandler(sys.stdout),
        ],
    )

def load_concat_samples(
    merged_data_paths: List[str], gtf_paths: List[str]
) -> Tuple[pd.DataFrame, List[str]]:
    samples = []
    samples_names = []
    for path, gtf_path in zip(merged_data_paths, gtf_paths):
        name = os.path.basename(path).split(".")[0]  # get name like sample_rep1
        samples_names.append(name)
        logging.info(f"Start processing {name}")
        samples.append(parse_sample(path, name, gtf_path))
    return (pd.concat(samples), samples_names)

def merge_count_file(
    count_files: List[str],
    label: str,
    transcript2gene: Dict[str, str],
    gene2symbol: Dict[str, str],
) -> pd.DataFrame:

    def parse_count_file(
        file_path: str,
        transcript2gene: Dict[str, str],
        gene2symbol: Dict[str, str],
        name: str,
    ) -> pd.DataFrame:
        count_df = pd.read_csv(file_path, sep="\t")
        count_df.columns = ["isoform_id", "count"]
        count_df["gene_id"] = count_df["isoform_id"].apply(
            lambda x: transcript2gene[x] if x in transcript2gene.keys() else np.nan
        )
        count_df["gene_symbol"] = count_df["gene_id"].apply(
            lambda x: gene2symbol[x] if x in gene2symbol.keys() else x
        )
        count_df[name] = count_df["count"]
        #     count_df = count_df.query('not feature_id.str.startswith("_")')
        logging.info(f"Loaded {name} count num")
        return count_df.query('not isoform_id.str.startswith("_")')

    count_dfs = []
    count_rep_names = []
    for count_file_path in count_files:
        rep_name = os.path.basename(count_file_path).split(".")[0]
        count_rep_names.append(rep_name)
        count_dfs.append(
            parse_count_file(count_file_path, transcript2gene, gene2symbol, rep_name)
        )
    merged_count_df = pd.DataFrame()
    merged_count_df = functools.reduce(
        lambda left, right: pd.merge(left, right, on=["isoform_id"], how="outer"),
        count_dfs,
    )
    merged_count_df.fillna(0, inplace=True)

    count_name = label + "_bam_count"
    merged_count_df[count_name] = merged_count_df[count_rep_names].sum(axis=1)
    cpm_name = label + "_bam_cpm"
    merged_count_df[cpm_name] = (
        merged_count_df[count_name] / merged_count_df[count_name].sum(axis=0) * 1000000
    )

    if len(count_files) > 1:
        return unify_columns(
            merged_count_df, "gene_id", "gene_symbol", ["_x", "_y"]
        ).iloc[
            :, [0, -4, -3, -2, -1]
        ]  # isoform_id, bam_count, bam_cpm, gene_id, gene_symbol
    else:
        return merged_count_df.iloc[:, [0, 2, 3, 5, 6]]


def generate_gtf_paths(paths: List[str]) -> List[str]:
    return [isoquant_derive_path(path, "gtf") for path in paths]


def generate_count_paths(paths: List[str]) -> List[str]:
    return [isoquant_derive_path(path, "count") for path in paths]


def isoquant_derive_path(path: str, file_type: str) -> str:
    assert file_type in ["gtf", "count"], "file_type should be gtf or count"
    base_folder_name = os.path.basename(os.path.normpath(path))
    if file_type == "gtf":
        return os.path.join(path, f"{base_folder_name}.extended_annotation.gtf")
    elif file_type == "count":
        return os.path.join(path, f"{base_folder_name}.transcript_counts.tsv")
    else:
        logging.info("invalid file_type, please check!")
    return None


def generate_geneview_df(readinfo_df: pd.DataFrame, label: str) -> pd.DataFrame:
    logging.info("Processing genes")
    sample_geneview = readinfo_df.loc[
        readinfo_df.polya_length > 10, ["polya_length", "isoform_id", "read_core_id"]
    ].set_index("isoform_id")
    sample_geneview = sample_geneview.groupby("isoform_id").agg(
        {
            "polya_length": [
                "count",
                "median",
                "mean",
                lambda x: x.quantile(0.25),
                lambda x: x.quantile(0.75),
            ]
        }
    )
    sample_geneview.columns = [
        label + "_fulllength_polya_read_num",
        label + "_fulllength_median_polya",
        label + "_fulllength_mean_polya",
        label + "_fulllength_25th",
        label + "_fulllength_75th",
    ]
    sample_geneview[f"{label}_fulllength_IQR"] = (
        sample_geneview[f"{label}_fulllength_75th"]
        - sample_geneview[f"{label}_fulllength_25th"]
    )
    return sample_geneview.reset_index()


def main():
    parser = argparse.ArgumentParser(description="Get gene cpm matrix for each samples")
    parser.add_argument(
        "--merged_data",
        nargs="+",
        required=True,
        help="Path to the merged Parquet file(s)",
    )
    parser.add_argument("--label", required=True, help="label for samples")
    parser.add_argument(
        "--output_dir", required=True, help="Directory to save analysis results"
    )
    parser.add_argument(
        "--isoquant_path", nargs="+", required=True, help="Path to isoquant dir"
    )

    args = parser.parse_args()

    merged_data_paths = args.merged_data
    output_dir = args.output_dir
    label = args.label
    isoquant_path = args.isoquant_path

    gtf_paths = generate_gtf_paths(isoquant_path)
    count_paths = generate_count_paths(isoquant_path)

    # Configure logging
    log_file_path = os.path.join(output_dir, f"polya_matrix_{label}.log")
    configure_logging(log_file_path)

    logging.info(f"Start processing {label} for cpm_gene_matrix")

    logging.info("Loading samples ...")
    rep_merged_df, samples_names = load_concat_samples(merged_data_paths, gtf_paths)

    # parse gtf
    logging.info(f"Processing gtf files ...")
    transcript_to_gene = {}
    gene_to_symbol = {}
    for gtf_path in gtf_paths:
        _, t2g, g2s = parse_gtf(gtf_path)
        transcript_to_gene.update(t2g)
        gene_to_symbol.update(g2s)

    # parse count file
    logging.info(f"Loading isoquant count files ...")
    merge_count = merge_count_file(
        count_paths, label, transcript_to_gene, gene_to_symbol
    )
    logging.info(f"merge_count cols: {merge_count.columns}")

    # process isoforms
    logging.info("Merging readinfo and count file")
    gene2polya = generate_geneview_df(rep_merged_df, label).merge(
        merge_count, on="isoform_id", how="outer"
    )
    gene2polya.fillna(0, inplace=True)

    gene2polya.to_csv(
        os.path.join(output_dir, f"{label}_cpm_polya_gene.csv"), sep="\t", index=None
    )
    logging.info(
        f"Successfully write output to {output_dir}/{label}_cpm_polya_gene.csv\n==============="
    )


if __name__ == "__main__":
    main()
