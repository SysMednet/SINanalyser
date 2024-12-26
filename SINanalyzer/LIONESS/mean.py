#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 16:17:24 2024

@author: tashatu
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 

def mean_calculation(data_df):
    print('Step 1 : mean calculation')
    if not isinstance(data_df, pd.DataFrame):
        raise ValueError("`data_df` must be a Pandas DataFrame with genes as rows and samples as columns.")

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
    for idx in range(data_sample_l): 
        edge_score = np.delete(data_val,idx,1)
        edge_score = np.corrcoef(edge_score)
        edge_score = edge_score[tri]
        edge_score[np.isnan(edge_score)] = 0
        edge_score =  data_sample_l * (agg - edge_score) + edge_score
        mean_sum += np.sum(edge_score)
        pair_num += pair_l
    mean = mean_sum/pair_num
    print('-> Mean of edge scores calculation ... Done')
    print(f'-> Mean = {mean}')
    return mean
