#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 21:58:23 2024

@author: tashatu
"""

import pandas as pd
import numpy as np

def correlation_calculation(data_df):
    print('Step 1 : weight calculation')
    data_val = data_df.values.T
    data_sample = list(data_df.columns)
    sample_corr = np.corrcoef(data_val)
    
    return sample_corr, data_sample

def weight_calculation(sample_corr,data_sample,outfile):
    sample_corr = (np.sum(sample_corr,axis=1)-1)/(len(data_sample)-1)
    corr_max , corr_min = np.max(sample_corr) , np.min(sample_corr)
    diff = corr_max - corr_min + 0.01
    weight = (sample_corr - corr_min + 0.01)/diff
    
    # Outfile
    sample_weight = {}
    with open(outfile, 'w') as file:
        file.write(f'Sample\tWeight\n')
        for s, w in zip(data_sample, weight):
            file.write(f'{s}\t{w}\n')
            sample_weight[s] =  w
    print('-> Sample weight calculation ... Done')
    return sample_weight
        
    