#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 02:19:25 2024

@author: tashatu
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 

def std_calculation(data_df,sample_weight,mean,amp=0.1):
    print('Step 3 : standard calculation')
    if not isinstance(data_df, pd.DataFrame):
        raise ValueError("`data_df` must be a Pandas DataFrame with genes as rows and samples as columns.")
    if not all(sample in data_df.columns for sample in sample_weight.keys()):
        raise ValueError("All keys in `sample_weight` must match the sample names in `data_df.columns`.")
    if not isinstance(amp, (int, float)):
        raise ValueError("`amp` must be a numeric value (int or float).")
    if not isinstance(mean, (int, float)):
        raise ValueError("`mean` must be a numeric value (int or float).")
        
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
    
    std_sum , pair_num = 0 , 0
    for idx,sample in enumerate(data_sample): 
        edge_score = np.corrcoef(np.c_[data_val,data_val[:,idx]])
        edge_score = edge_score[tri]
        edge_score[np.isnan(edge_score)] = 0
        edge_score =  sample_weight[sample] * amp * data_sample_l * (edge_score - agg) + agg
        edge_score = (edge_score-mean)**2
        std_sum += np.sum(edge_score)
        pair_num += pair_l
    std = np.sqrt(std_sum/(pair_num-1))
    print('-> Standard of edge scores calculation ... Done')
    print(f'-> Standard = {std}')
    return std
    