#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Dict, List
import pandas as pd
import numpy as np
import os
import sys
import gzip

def parse_gtf(file_path):
    data = {
        'seqname': [],
        'source': [],
        'feature': [],
        'start': [],
        'end': [],
        'score': [],
        'strand': [],
        'frame': [],
        'attribute': [],
        'gene_id': [],
        'transcript_id': [],
        'gene_name': [],
        'mgi_id': []
    }

    transcript_to_gene = {}
    gene_to_symbol = {}

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue 

            attributes = fields[8]

            attribute_dict = {}
            for attr in attributes.split(';'):
                if attr.strip():
                    key_value = attr.strip().split(' ')
                    if len(key_value) >= 2:
                        key = key_value[0]
                        value = ' '.join(key_value[1:]).strip('"')
                        attribute_dict[key] = value

            gene_id = attribute_dict.get('gene_id', '')
            transcript_id = attribute_dict.get('transcript_id', '')
            gene_name = attribute_dict.get('gene_name', '')
            mgi_id = attribute_dict.get('mgi_id', '')

            data['seqname'].append(fields[0])
            data['source'].append(fields[1])
            data['feature'].append(fields[2])
            data['start'].append(fields[3])
            data['end'].append(fields[4])
            data['score'].append(fields[5])
            data['strand'].append(fields[6])
            data['frame'].append(fields[7])
            data['attribute'].append(attributes)
            data['gene_id'].append(gene_id)
            data['transcript_id'].append(transcript_id)
            data['gene_name'].append(gene_name)
            data['mgi_id'].append(mgi_id)

            if transcript_id:
                transcript_to_gene[transcript_id] = gene_id
            if gene_name:
                gene_to_symbol[gene_id] = gene_name

    df = pd.DataFrame(data)

    return df, transcript_to_gene, gene_to_symbol

def read_read_info(file_path):
    try:
        return pd.read_csv(
            file_path, 
            sep='\t', 
            # chunksize=chunksize, 
            dtype={"retention_introns": str},
            # low_memory=False 
        )
    except Exception as e:
        print(f"Error reading read info file {file_path}: {e}", file=sys.stderr)
        return iter([])

def optimize_dataframe(df):
    for col in df.columns:
        col_data = df[col]
        if col_data.dtype == 'object':
            try:
                col_data_numeric = pd.to_numeric(col_data, errors='ignore')
                if col_data_numeric.dtype != 'object':
                    df[col] = col_data_numeric
                else:
                    df[col] = df[col].astype(str)
            except Exception as e:
                df[col] = df[col].astype(str)
        elif pd.api.types.is_integer_dtype(col_data):
            df[col] = pd.to_numeric(col_data, downcast='integer')
        elif pd.api.types.is_float_dtype(col_data):
            df[col] = pd.to_numeric(col_data, downcast='float')

    for col in df.select_dtypes(include=['object']).columns:
        num_unique_values = df[col].nunique()
        num_total_values = len(df[col])
        if num_unique_values / num_total_values < 0.5:
            df[col] = df[col].astype('category')

    return df

def merge_dict(d1: Dict[str, str], d2: Dict[str, str]) -> Dict[str, str]:
    merged_dict = {**d1, **d2}
    return {key: merged_dict[key] for key in set(merged_dict)}


def unify_columns(df: pd.DataFrame, colA: str, colB: str, suffix: List[str] = ['_x', '_y']) -> pd.DataFrame:
    '''
    description: due with col_x and col_y after pd.merge
    return dataframe
    '''
    def merge_columns_vectorized(cols: List[str]) -> pd.DataFrame:
        col_x, col_y = cols
        same_values = df[col_x] == df[col_y]
        different_values = ~same_values & ~((df[col_x] == 0) | (df[col_x].isna()) | (df[col_y] == 0) | (df[col_y].isna()))

        if different_values.any():
            raise Exception(f"Bug detected: {col_x} and {col_y} have different non-zero values.")

        merged_col = df[col_x].where(same_values | (df[col_x] != 0) & df[col_x].notna(), df[col_y])
        return merged_col

    df[colA] = merge_columns_vectorized([colA+s for s in suffix])
    df[colB] = merge_columns_vectorized([colB+s for s in suffix])

    return df