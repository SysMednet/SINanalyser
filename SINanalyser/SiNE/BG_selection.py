#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 17:42:45 2024

@author: tashatu
"""

import sys
import pandas as pd
import random

random.seed(42)

def bin_choose_nonselect(df, select_gene, n_genes, fold):
    """
    Helper function to choose non-selected genes for background based on binning strategy.

    Parameters:
    df (pandas.DataFrame): The input gene expression DataFrame.
    select_gene (list): List of selected gene names.
    n_genes (int): Number of selected genes.
    fold (int): Fold parameter for selecting non-selected genes.

    Returns:
    list: Indices of the chosen non-selected genes.
    """
    bin_dist = []
    gene_indices = []
    # Separate bins and calculate frequency
    mean = df.mean(axis = 1, skipna = False)
    sort = pd.Series(mean.sort_values(ascending=False))
    bins = pd.cut(sort, bins=10, labels=False, duplicates='drop')
    df_bin_count = pd.DataFrame(bins.value_counts().sort_index())
    df_bin_count.columns = ['Counts']
    df_bin_count['Frequency'] = round(df_bin_count['Counts']/len(df.index),3)
    df_mean_bin = pd.DataFrame(bins, columns=['Bin'])
    df_mean_bin['Select'] = df_mean_bin.index.isin(select_gene).astype(int)
    total_genes = n_genes*(fold+1)
    # Pick n_genes*fold*frequency non-selected genes from each bin
    for idx, row in df_bin_count.iterrows():
        # Calculate the number of non-selected genes should be chosen
        bin_gene_num = round(total_genes * row['Frequency']) - df_mean_bin[(df_mean_bin['Bin'] == idx) & (df_mean_bin['Select'] == 1)].shape[0]
        # Random choose from non-selected genes
        if bin_gene_num > 0: # if <= 0 , should not choose any genes
            bin_samples = list(df_mean_bin[(df_mean_bin['Bin'] == idx) & (df_mean_bin['Select'] == 0)].index)
            if len(bin_samples) < bin_gene_num: bin_gene_num = len(bin_samples)
            bin_samples_choose = random.sample(bin_samples, k=bin_gene_num)
            gene_indices.extend(bin_samples_choose)
    return gene_indices


def background_GEM_filter(df, filter_GEM, fold=5, outdir=None):
    """
    Filters the gene expression matrix (GEM) to include selected genes and a background of non-selected genes.
    
    Parameters:
    df (pandas.DataFrame): The input gene expression DataFrame. The data should be formatted with genes as rows and samples as columns.
    filter_GEM (pandas.DataFrame): DataFrame of selected genes.
    fold (int): The fold parameter for the number of non-selected genes to include as background (default=5).
    outdir (str or None): Output directory to save the background GEM. If None, no file will be saved.
    
    Returns:
    pandas.DataFrame: Background gene expression matrix containing both selected and non-selected genes.
    
    Parameter
    ----------
    df : pandas.DataFrame
        The input gene expression DataFrame with genes as rows and samples as columns.
    
    filter_GEM : pandas.DataFrame
        DataFrame of selected genes to include in the output.
    
    fold : int, optional
        The fold parameter to control the ratio of non-selected genes to selected genes (default=5).
    
    outdir : str or None, optional
        Output directory to save the background GEM file. If None, the file will not be saved.
    
    Notes
    -----
    The function ensures no duplicate genes by averaging their values if they appear multiple times.
    Selected genes are retained in the background GEM, and non-selected genes are chosen based on bin frequency to match the desired ratio.
    """
    df = df.groupby(level=0).mean() # to avoid duplicate genes
    select_gene = list(filter_GEM.index) # number of selected-genes
    n_genes = len(select_gene) # number of selected-genes
    nonselect_gene = bin_choose_nonselect(df, select_gene, n_genes, fold)
    gene_total = list(select_gene + nonselect_gene)# select, non-select
    BG_GEM = df.loc[gene_total]
    if outdir: BG_GEM.to_csv(f'{outdir}/BG_data.txt', sep="\t", header=True, index=True)
    return BG_GEM