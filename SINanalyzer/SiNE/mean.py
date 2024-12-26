#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 22:59:56 2024

@author: tashatu
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 

def mean_calculation(data_df,sample_weight,amp=0.1):
    print('Step 2 : mean calculation')
    if not isinstance(data_df, pd.DataFrame):
        raise ValueError("`data_df` must be a Pandas DataFrame with genes as rows and samples as columns.")
    if not all(sample in data_df.columns for sample in sample_weight.keys()):
        raise ValueError("All keys in `sample_weight` must match the sample names in `data_df.columns`.")
    if not isinstance(amp, (int, float)):
        raise ValueError("`amp` must be a numeric value (int or float).")
        
    # Data
    data_gene = list(data_df.index)
    data_sample = list(data_df.columns)
    data_val = data_df.values
    data_gene_l = len(data_gene)
    data_sample_l = len(data_sample)
    print('-> Data preprocessing ... Done')
    
    # Aggregate network construction
    agg = np.corrcoef(data_val)
    tri = np.full(agg.shape,False)
    tri[np.triu_indices(data_gene_l,1)] = True
    agg = agg[tri]
    agg[np.isnan(agg)] = 0
    pair_l = len(agg)
    print('-> Aggregate network construction ... Done')
    
    mean_sum , pair_num = 0 , 0
    for idx,sample in enumerate(data_sample): 
        edge_score = np.corrcoef(np.c_[data_val,data_val[:,idx]])
        edge_score = edge_score[tri]
        edge_score[np.isnan(edge_score)] = 0
        edge_score =  sample_weight[sample] * amp * data_sample_l * (edge_score - agg) + agg
        mean_sum += np.sum(edge_score)
        pair_num += pair_l
    mean = mean_sum/pair_num
    print('-> Mean of edge scores calculation ... Done')
    print(f'-> Mean = {mean}')
    return mean
