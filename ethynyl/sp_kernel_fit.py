#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HCGA based graph

- Builds weighted graphs 
- puts betweenness centrality as nodes 
- computes the shortest path kernel 
"""
#%%
import numpy as np
import pandas as pd
import os, glob, re, sys
import networkx as nx
from scipy.spatial import distance_matrix
import subprocess as sp
from sklearn.preprocessing import MinMaxScaler

import grakel 
import pickle 
from grakel import Graph
from grakel.kernels import PyramidMatch, ShortestPathAttr, WeisfeilerLehman, ShortestPath
import time
import os
import numpy as np
import csv

"""       BUILDS A BINARY OR WEIGHTED GRAPH FROM XYZ FILE


        Possible Inputs: 3x3x3 Supercell, NeighShell (15 molecules)

        Graph Edges seen as contacts distance cutoff: sum(vdw_radii)+2 A
        Edges = Distances ≤ sum of the VdW radii + 2 Angstrom will be considered
        (CCDC_Compack inspired cutoff)

"""

# ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹ USER INPUT ››››››››››››››››››››››››››››››››››››››››››››››››››
inp_folder = 'Neighshells'#Neighshells'   # dir with xyz to convert into graphs


## variables which can be easily changed, hard coded for comparison of graphs


atom_vdw_radii = {
              'Al': 2, 'Sb': 2, 'Ar': 1.88, 'As': 1.85, 'Ba': 2,
              'Be': 2, 'Bi': 2, 'B': 2, 'Br': 1.85, 'Cd': 1.58,
              'Cs': 2, 'Ca': 2, 'C': 1.7, 'Ce': 2, 'Cl': 1.75,
              'Cr': 2, 'Co': 2, 'Cu': 1.4, 'Dy': 2, 'Er': 2,
              'Eu': 2, 'F':  1.47, 'Gd': 2, 'Ga': 1.87, 'Ge': 2,
              'Au': 1.66, 'Hf': 2, 'He': 1.4, 'Ho': 2, 'H': 1.09,
              'In': 1.93, 'I': 1.98, 'Ir': 2, 'Fe': 2, 'Kr': 2.02,
              'La': 2, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2, 'Mg': 1.73,
              'Mn': 2, 'Hg': 1.55, 'Mo': 2, 'Nd': 2, 'Ne': 1.54,
              'Ni': 1.63, 'Nb': 2, 'N':  1.55, 'Os': 2, 'O':  1.52,
              'Pd': 1.63, 'P': 1.8, 'Pt': 1.72, 'K': 2.75, 'Pr': 2,
              'Pa': 2, 'Re': 2, 'Rh': 2, 'Rb': 2, 'Ru': 2, 'Sm': 2,
              'Sc': 2, 'Se': 1.9, 'Si': 2.1, 'Ag': 1.72, 'Na': 2.27,
              'Sr': 2, 'S': 1.8, 'Ta': 2, 'Te': 2.06, 'Tb': 2,
              'Tl': 1.96, 'Th': 2, 'Tm': 2, 'Sn': 2.17, 'Ti': 2,
              'W': 2, 'U':  1.86, 'V':  2, 'Xe': 2.16, 'Yb': 2,
              'Y': 2, 'Zn': 1.29, 'Zr': 2, 'X':  1.0, 'D':  1.0
                 }



skip = 0 #2


def build_graph(inp_folder, basename, labels, coords, atom_vdw_radii, weighted=True, sparsify=False):
    # turns distance matrix into a simple 0 1 binary graph according to
    # sum of VdW distance + X A. see header
    # get vdw radii for each atom to compute sum of vdw radii
    dist_mat = distance_matrix(coords, coords)
    add_to_threshold = 2  # extra 2 A
    # skip = 2

    labels_vdw= []
    for i in range(0,len(labels)):
        vdw = atom_vdw_radii.get(str(labels[i]))
        labels_vdw.append(vdw)
    # turning distance matrix into an unweighted graph - simplest version
    # decide what kind of distance matrix you want:
    counter=0

    if weighted==False:
        """ BINARY GRAPH """
        for i in range(0,int(np.shape(labels)[0])):
            for j in range(0,int(np.shape(labels)[0])):
                sum_vdw_radii = labels_vdw[i]+labels_vdw[j]
                threshold = float(sum_vdw_radii) +2

                if float(dist_mat[i, j]) <= threshold:
                        dist_mat[i,j]= 1
                else:
                    dist_mat[i,j]= 0

                counter =counter + 1


    else:
        """ WEIGHTED GRAPH """
        for i in range(0,int(np.shape(labels)[0])):
            for j in range(0,int(np.shape(labels)[0])):
                sum_vdw_radii = labels_vdw[i]+labels_vdw[j]
                threshold = float(sum_vdw_radii) + add_to_threshold

                if float(dist_mat[i, j]) >= threshold:
                    dist_mat[i,j]= 0

                counter =counter + 1


    
    
    feats = pd.get_dummies(pd.DataFrame(labels)).values
    

    G = nx.from_numpy_matrix(dist_mat)
    for i in range(len(G)):
        G.nodes[i]['atom_types'] = feats[i,:]


    

    return G


#%%
#‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹››››››››››››››››››››››››››››››››››››››


def construct_helicene_graphs():
    
    skip = 2
    inp_folder = 'Neighshells'#Neighshells'   # dir with xyz to convert into graphs
    # create output folder if it doesn't already exist
    # if os.path.isdir(f"csv_{inp_folder}") == True:
    #     pass
    # else: 
    #     sp.call(f"mkdir csv_{inp_folder}", shell=True)
    
    
    open("HELICENES_A.txt", "w+") # contains all the adjacency matrices/graphs (X)
    open("HELICENES_graph_indicator.txt", "w+")
    open("HELICENES_graph_labels.txt", "w+") # contains the y target var
    open("HELICENES_node_attributes.txt", "w+") # contains the atom labels
    count = 0
    
    graphs = []
    graph_labels = []
    node_attr = []
    names_list = []
    grakel_graphs = []
    # loop through input folder
    for name in glob.glob(f"./{inp_folder}/*.xyz"):
        print(name)
        basename = re.findall("\d+", name)[-1] #last item get right number if res3_supercell_0.xyz => extracts 0
        print(basename)
        count = count+1
        
        
        ids = re.findall("\d+", name) 
    
        names_list.append(basename)
        
        """ READ INPUT FILES """
        labels = np.genfromtxt(name,usecols=0,dtype=str, skip_header=skip)
        coords = np.genfromtxt(name, skip_header=skip)[:,1:]
    
        """ BUILD BINARY DISTANCE MATRIX """
        graph = build_graph(inp_folder, basename, labels, coords,atom_vdw_radii)
        
        G = edge_betweenness(graph)
        graphs.append(G)
        # print(len(graph))
        # node_attr.append(labels)
        # grakel_graphs.append(grakel_graph)
    return graphs, names_list #,graph_labels


# %%

def edge_betweenness(G):
    node_labels = nx.betweenness_centrality(G, normalized=True)
    nx.set_node_attributes(G, node_labels, 'label')
    return G

def convert_to_grakel_graphs(Gs):
    return grakel.graph_from_networkx(Gs, node_labels_tag='label')

# def read_in_grakel(graph, node_attr):
#     # read graph from pickle
#     G = graph#nx.from_numpy_matrix(graph)#nx.read_gpickle(f"{working_dir}/{folder}/{filn}.pkl")
#     # compute and add node_labels on betweenness centrality
#     # node_labels = node_attr#nx.betweenness_centrality(G, normalized=True)
#     node_labels = nx.betweenness_centrality(G, normalized=True)
#     nx.set_node_attributes(G, node_labels, 'label')
#     # turn nx graph into grakel graph
#     return grakel.graph_from_networkx([G], node_labels_tag='label')


def real_run(kernel_name, kernel, graphs, names_list):
    pm = kernel

    start2 = time.time()
    # sim_mat = np.empty((n_dim, n_dim))
    matrix = []
    graphs_list = list(graphs)
    similarity = pm.fit_transform(graphs_list)
    # print(type(similarity))
    # print(similarity)

    np.savetxt(f"{kernel_name}.csv", similarity, delimiter=",")

    # with open(f"{kernel_name}.csv", "w", newline="") as f:
    #     writer = csv.writer(f)
    #     writer.writerows(matrix)

    end2 = time.time()
    print("Time needed for Kernel comparison (s): ", end2 - start2)


# make graphs from xyz
graphs, names_list = construct_helicene_graphs()
grakel_graphs = convert_to_grakel_graphs(graphs)

sp = ShortestPath()
real_run(f"SP", sp, grakel_graphs, names_list)
