#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 18:04:24 2024

@author: tashatu
"""
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 
from .SiNE.BG_selection import bin_choose_nonselect, background_GEM_filter
from .SiNE.weight import correlation_calculation as SiNE_correlation, weight_calculation as SiNE_weight
from .SiNE.mean import mean_calculation as SiNE_mean
from .SiNE.std import std_calculation as SiNE_std
from .SiNE.network import network_construction as SiNE_network
from .SWEET.weight import correlation_calculation as SWEET_correlation, weight_calculation as SWEET_weight
from .SWEET.mean import mean_calculation as SWEET_mean
from .SWEET.std import std_calculation as SWEET_std
from .SWEET.network import network_construction as SWEET_network
from .LIONESS.mean import mean_calculation as LIONESS_mean
from .LIONESS.std import std_calculation as LIONESS_std
from .LIONESS.network import network_construction as LIONESS_network
from .SSN.network import network_construction as SSN_network


def log2_normalization(df):
    """
    Applies log2 normalization to the gene expression DataFrame.

    Parameters:
    df (pandas.DataFrame): The input gene expression DataFrame, with genes as rows and samples as columns.

    Returns:
    pandas.DataFrame: Log2-normalized gene expression DataFrame..
    """
    df_log2 = df + 1
    df_log2 = np.log2(df_log2)
    return df_log2

def data_preprocessing(data, gene_set=None, normalize=False):
    """
    Preprocess the input gene expression data.

    Parameters:
    data (pd.DataFrame): The gene expression data.
    gene_set (list or None): A list of genes to focus on.
    normalize (bool): Whether to apply log2 normalization.

    Returns:
    pd.DataFrame: The preprocessed data.
    """
    if normalize:
        data = log2_normalization(data)

    if gene_set is not None:
        filter_data = data.loc[data.index.intersection(gene_set)]
    else:
        filter_data = None
    return data, filter_data

def network_construction(case_data, control_data=None, method="SiNE", amp=0.1, pvalue=0.05, gene_set=None, sample_list=None, normalize=False, outdir=None, output_format="npz"):
    """
    Executes the network construction using the specified method (SiNE, LIONESS, or SSN).

    Parameters:
    case_data (pd.DataFrame): A DataFrame containing the case gene expression data. The data should be formatted with genes as rows and samples as columns.
    control_data (pd.DataFrame or None): A DataFrame containing the control gene expression data (only required for SSN). If not used, default is None.
    method (str): The network construction method to use. Options are "SiNE", "SWEET", "LIONESS", or "SSN".
    amp (float): Amplitude parameter for SiNE and SWEET network construction (default=0.1).
    pvalue (float): P-value threshold for network construction (default=0.05).
    gene_set (list or None): A list of genes to focus on. If None, all genes in the dataset are used.
    sample_list (list or None): A list of sample names to include in the output. If None, all samples are output.
    normalize (bool): If True, the data will be log2 normalized (default=False).
    outdir (str or None): Output directory to save the results. If None, the current directory is used.
    output_format (str): Output file format. Options are "npz", "edge_list", "edge_list_score", "edge_list_zscore".

    Returns:
    None

    Notes:
    This function dynamically adjusts parameters and execution based on the specified method.
    """
    if outdir is None:
        outdir = os.getcwd()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Load data
    if not isinstance(case_data, pd.DataFrame):
        raise ValueError("case_data must be a Pandas DataFrame.")
    case_data, filter_data = data_preprocessing(case_data, gene_set, normalize)

    # SiNE
    if method == "SiNE":
        if control_data is not None:
            raise ValueError("Data control should not be provided when method is 'SiNE'.")
        ################################ Sample weight ################################
        sample_corr, data_sample = SiNE_correlation(case_data)
        weight_file = f'{outdir}/weight.txt'
        sample_weight = SiNE_weight(sample_corr,data_sample,weight_file)
        if gene_set != None:
            ################################ Background genes ################################
            if not filter_data.empty:
                BG_data = background_GEM_filter(case_data, filter_data, outdir=outdir)
            ################################ Mean, Standard, Network ################################
            mean = SiNE_mean(BG_data,sample_weight)
            std = SiNE_std(BG_data,sample_weight,mean)
            SiNE_network(filter_data,sample_weight,mean,std,amp=amp,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
        else:
            ################################ Mean, Standard, Network ################################
            mean = SiNE_mean(case_data,sample_weight)
            std = SiNE_std(case_data,sample_weight,mean)
            SiNE_network(case_data,sample_weight,mean,std,amp=amp,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
            
    # SWEET
    elif method == "SWEET":
        if control_data is not None:
            raise ValueError("Data control should not be provided when method is 'SWEET'.")
        ################################ Sample weight ################################
        sample_corr, data_sample = SWEET_correlation(case_data)
        weight_file = f'{outdir}/weight.txt'
        sample_weight = SWEET_weight(sample_corr,data_sample,weight_file)
        if gene_set != None:
            ################################ Mean, Standard, Network ################################
            mean = SWEET_mean(filter_data,sample_weight)
            std = SWEET_std(filter_data,sample_weight,mean)
            SWEET_network(filter_data,sample_weight,mean,std,amp=amp,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
        else:
            ################################ Mean, Standard, Network ################################
            mean = SWEET_mean(case_data,sample_weight)
            std = SWEET_std(case_data,sample_weight,mean)
            SWEET_network(case_data,sample_weight,mean,std,amp=amp,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
    
    # LIONESS
    elif method == "LIONESS":
        if control_data is not None:
            raise ValueError("Data control should not be provided when method is 'LIONESS'.")
        if gene_set != None:
            ################################ Mean, Standard, Network ################################
            mean = LIONESS_mean(filter_data)
            std = LIONESS_std(filter_data,mean)
            LIONESS_network(filter_data,mean,std,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
        else:
            ################################ Mean, Standard, Network ################################
            mean = LIONESS_mean(case_data)
            std = LIONESS_std(case_data,mean)
            LIONESS_network(case_data,mean,std,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
        
    # SSN 
    elif method == "SSN":
        if control_data is None:
            raise ValueError("Control dataset must be provided for SSN method.")
        if not isinstance(control_data, pd.DataFrame):
            raise ValueError("control_data must be a Pandas DataFrame.")
        control_data, control_filter_data = data_preprocessing(control_data, gene_set, normalize)
        if gene_set != None:
            ################################ Network ################################    
            SSN_network(filter_data,control_filter_data,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
        else:
            ################################ Network ################################    
            SSN_network(case_data,control_data,pvalue=pvalue,sample_list=sample_list,outdir=outdir,output_format=output_format)
    else:
        raise ValueError("The network construction method must be 'SiNE', 'SWEET', 'LIONESS', or 'SSN'")
    
