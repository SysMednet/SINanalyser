#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 15:48:51 2024

@author: tashatu
"""
from pyvis.network import Network
import scipy.sparse
from scipy.sparse import coo_matrix

# Step 1 (option 1): Read edge_list.npz 
def read_npz(file_path, gene_label):
    """
    Reads a .npz file containing the adjacency matrix and returns edges with gene labels.

    Parameters:
    file_path (str): Path to the .npz file containing the adjacency matrix.
    gene_label (str): Path to the gene label file mapping indices to gene names.

    Returns:
    list: A list of tuples representing edges in the format (source, target).

    Example:
    edges = read_npz("matrix.npz", "gene_label.txt")
    """
    genes = {}
    edges = []
    with open(gene_label,mode='r') as file :
        for line in file :
            idx, name = line.strip('\n').split('\t')
            genes[idx] = name
    matrix = scipy.sparse.load_npz(file_path).toarray()
    load_sparse = coo_matrix(matrix)
    for i,j,k in zip(load_sparse.row,load_sparse.col,load_sparse.data) :
        edges.append((genes[str(i)], genes[str(j)]))
    return edges

# Step 1 (option 2): Read edge_list.txt 
def read_edge_list(file_path):
    """
    Reads an edge list file and returns edges.

    Parameters:
    file_path (str): Path to the edge list file.

    Returns:
    list: A list of tuples representing edges in the format (source, target).

    Example:
    edges = read_edge_list("edge_list.txt")
    """
    edges = []
    with open(file_path, 'r') as file:
        for line in file:
            nodes = line.strip().split('\t')
            if len(nodes) == 2:  # two nodes (source, target)
                edges.append((nodes[0], nodes[1]))
    return edges

# Step 1 (option 3): Read edge_list_score.txt or  edge_list_zscore.txt
def read_edge_list_with_score(file_path):
    """
    Reads an edge list file with scores and returns edges with their weights.

    Parameters:
    file_path (str): Path to the edge list file.

    Returns:
    list: A list of tuples representing edges in the format (source, target, score).

    Example:
    edges = read_edge_list_with_score("edge_list.txt")
    """
    edges = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) == 3:  # source, target, score
                source, target, score = parts
                edges.append((source, target, float(score)))
    return edges

# Step 2: Using Pyvis to create the interective network
def create_pyvis_network(edges, output_html):
    """
    Creates an interactive network visualization using Pyvis and saves it to an HTML file.

    Parameters:
    edges (list): A list of tuples representing edges in the format (source, target).
    output_html (str): Path to save the HTML file.

    Returns:
    None

    Example:
    create_pyvis_network(edges, "network.html")
    """
    net = Network(height='600px', width='80%', bgcolor='#ffffff', # Changed height
                font_color='black',notebook = True, directed=False, select_menu=True, cdn_resources='in_line')
    
    for edge in edges:
        source, target, score = edge
        net.add_node(source, label=source)
        net.add_node(target, label=target)
        if score == None:
            net.add_edge(source, target)
        else:
            net.add_edge(source, target, value=score)  # Use 'value' for edge thickness

    net.toggle_physics(False)
    net.show_buttons(filter_=['physics'])
    net.save_graph(output_html)
    print(f"The network graph has been saved to {output_html}")

def graph_plot(file_path, gene_label=None, save_file_path="network.html", input_format="npz"):
    """
    Visualizes a network from either a .npz file or an edge list file and saves it as an HTML file.

    Parameters:
    file_path (str): Path to the input file (.npz or edge list).
    gene_label (str or None): Path to the gene label file (required for .npz format).
    save_file_path (str): Path to save the output HTML file.
    input_format (str): Input file format. Options are "npz" or "edge_list".

    Returns:
    None

    Examples:
    - From an .npz file:
      graph_plot(file_path="matrix.npz", gene_label="gene_label.txt", save_file_path="network.html", input_format="npz")

    - From an edge list file:
      graph_plot(file_path="edge_list.txt", save_file_path="network.html", input_format="edge_list")
    """
    # npz
    if input_format == "npz":
        if gene_label == None:
            raise ValueError("gene_label file should be uploaded.")
        edges = read_npz(file_path, gene_label)
        create_pyvis_network([(e[0], e[1], None) for e in edges], output_html=save_file_path)  
    # edge list
    elif input_format == "edge_list":
        if gene_label != None:
            raise ValueError("gene_label file should not be uploaded.")
        edges = read_edge_list(file_path)
        create_pyvis_network([(e[0], e[1], None) for e in edges], output_html=save_file_path)
    elif input_format == "edge_list_score" or input_format == "edge_list_zscore":
        if gene_label != None:
            raise ValueError("gene_label file should not be uploaded.")
        edges = read_edge_list_with_score(file_path)
        create_pyvis_network(edges, output_html=save_file_path)
    else:
        raise ValueError("Invalid input_format. Choose either 'npz', 'edge_list', 'edge_list_score' or 'edge_list_zscore'.")
        


