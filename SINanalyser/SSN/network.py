#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 17:01:29 2024

@author: tashatu
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 

def network_construction(data_df,data_control_df,pvalue=0.05,sample_list=None,outdir=None,output_format="npz"):
    if not isinstance(data_df, pd.DataFrame):
        raise ValueError("`data_df` must be a Pandas DataFrame with genes as rows and samples as columns.")
    if not isinstance(data_control_df, pd.DataFrame):
        raise ValueError("`data_control_df` must be a Pandas DataFrame with genes as rows and samples as columns.")
    if output_format not in ["npz", "edge_list", "edge_list_score", "edge_list_zscore"]:
        raise ValueError("`output_format` must be one of 'npz', 'edge_list', 'edge_list_score', or 'edge_list_zscore'.")
    if sample_list is not None and not all(sample in data_df.columns for sample in sample_list):
        raise ValueError("All samples in `sample_list` must exist in `data_df.columns`.")
        
    if not os.path.exists(outdir): os.makedirs(outdir) # Check if the directory exists

    zscore = stats.norm.isf(pvalue/2)
    # Data
    data_gene = list(data_df.index)
    data_sample = list(data_df.columns)
    data_val = data_df.values
    data_gene_l = len(data_gene)
    data_sample_l = len(data_sample)
    if sample_list == None: sample_list = data_sample
    # Data_control
    data_control_gene = list(set(data_control_df.index)&set(data_gene)) # overlap with data_gene
    data_control_sample = list(data_control_df.columns)
    data_control_df = data_control_df.loc[data_control_gene] # overlap with data_gene
    data_control_val = data_control_df.values
    data_control_sample_l = len(data_control_sample)
    if set(data_gene) != set(data_control_gene) :
        print('The genes in the two matrices are different.')
        return False
    print('-> Data preprocessing ... Done')
    
    if output_format == "npz":
        # Gene index mapping table
        with open(f'{outdir}/Gene_index_mapping_table.txt','w') as gene_file :
            for idx,val in enumerate(data_gene) :
                gene_file.write(f'{idx}\t{val}\n')
        print('-> Gene index mapping table ... Done')
    
    # Aggregate network construction with data_control
    agg = np.round(np.corrcoef(data_control_val),6)
    agg[agg==1] = 0.99999
    agg[agg==-1] = -0.99999
    agg[np.isnan(agg)] = 0
    agg_zero = np.where(agg == 0)
    agg_deno = (1-np.square(agg))/(data_control_sample_l-1)    
    agg_deno[agg_zero] = 1
    print('-> Aggregate network construction ... Done')
    
    # To store network in triangular matrix
    edge = np.full((data_gene_l,data_gene_l),False)
    tri = np.full((data_gene_l,data_gene_l),False)
    tri[np.triu_indices(data_gene_l,1)] = True

    # Network contruction
    for idx,sample in enumerate(data_sample):
        if sample in sample_list:
            edge_score = np.c_[data_control_val,data_val[:,idx]]
            edge_score = np.round(np.corrcoef(edge_score),6)
            edge_score[np.isnan(edge_score)] = 0
            edge_score_z = np.abs((edge_score-agg)/agg_deno)
            edge_score_z[range(data_gene_l), range(data_gene_l)] = 0 # Self correlation=0
            edge_score_z[agg_zero] = 0
            edge[tri] = edge_score_z[tri] >= zscore
            rows, cols = np.where(edge)
            
            # Handle output formats
            if output_format == "npz":
                sparse_s = csr_matrix(edge)
                scipy.sparse.save_npz(f'{outdir}/{sample}.npz', sparse_s)
     
            elif output_format == "edge_list":
                with open(f'{outdir}/{sample}_edge_list.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\n')

            elif output_format == "edge_list_score":
                with open(f'{outdir}/{sample}_edge_list_score.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\t{edge_score[i, j]:.6f}\n')

            elif output_format == "edge_list_zscore":
                with open(f'{outdir}/{sample}_edge_list_zscore.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\t{edge_score_z[i, j]:.6f}\n')
    print('-> Single-sample network construction ... Done')