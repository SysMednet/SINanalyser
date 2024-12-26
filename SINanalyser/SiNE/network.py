#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 02:32:46 2024

@author: tashatu
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 

def network_construction(data_df,sample_weight,mean,std,amp=0.1,pvalue=0.05,sample_list=None,outdir=None,output_format="npz"):
    print('Step 4 : network construction')
    if not isinstance(data_df, pd.DataFrame):
        raise ValueError("`data_df` must be a Pandas DataFrame with genes as rows and samples as columns.")
    if not all(sample in data_df.columns for sample in sample_weight.keys()):
        raise ValueError("All keys in `sample_weight` must match the sample names in `data_df.columns`.")
    if not isinstance(amp, (int, float)):
        raise ValueError("`amp` must be a numeric value (int or float).")
    if not isinstance(mean, (int, float)):
        raise ValueError("`mean` must be a numeric value (int or float).")
    if not isinstance(std, (int, float)):
        raise ValueError("`standard` must be a numeric value (int or float).")
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
    print('-> Data preprocessing ... Done')

    if output_format == "npz":
        # Gene index mapping table
        with open(f'{outdir}/Gene_index_mapping_table.txt','w') as gene_file :
            for idx,val in enumerate(data_gene) :
                gene_file.write(f'{idx}\t{val}\n')
        print('-> Gene index mapping table ... Done')
    
    # Aggregate network construction
    agg = np.corrcoef(data_val)
    tri = np.full(agg.shape,False)
    tri[np.triu_indices(data_gene_l,1)] = True
    agg = agg[tri]
    agg[np.isnan(agg)] = 0
    print('-> Aggregate network construction ... Done')
    
    # Network contruction
    edge = np.full((data_gene_l,data_gene_l),False)
    for idx,sample in enumerate(data_sample):
        if sample in sample_list:
            edge_score = np.corrcoef(np.c_[data_val,data_val[:,idx]])
            edge_score = edge_score[tri]
            edge_score[np.isnan(edge_score)] = 0
            edge_score =  sample_weight[sample] * amp * data_sample_l * (edge_score - agg) + agg
            edge_score_z = np.abs((edge_score-mean)/std)
            edge.fill(False)
            edge[tri] = edge_score_z >= zscore
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
                        score_index = data_gene_l*i-(i*(i+1))//2+(j-i-1)
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\t{edge_score[score_index]:.6f}\n')

            elif output_format == "edge_list_zscore":
                with open(f'{outdir}/{sample}_edge_list_zscore.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        score_index = data_gene_l*i-(i*(i+1))//2+(j-i-1)
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\t{edge_score_z[score_index]:.6f}\n')

    print('-> Single-sample network construction ... Done')
