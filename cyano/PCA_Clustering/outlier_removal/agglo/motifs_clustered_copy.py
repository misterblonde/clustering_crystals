#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

ANALYSIS OF GRAPH SIMILARITY KERNELS

1. PCA plot
2. Clustering
@author: JAS

"""

#%%
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns 
from adjustText import adjust_text
import sys
import sklearn
from sklearn.cluster import DBSCAN
from sklearn.datasets import load_iris
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
from statistics import mean 
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import v_measure_score
from sklearn.decomposition import PCA
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, dendrogram


def see_five_highest_sim(adj):
    """
    Function to sort top similarity scores in descending order of the entire 
    array, such that the number of duplicate packings present in the landscape 
    can be determined for rescaling. 

    """
    x = adj.values
    flat = x.flatten()
    flat = list(set(flat))
    flat.sort(reverse=True)
    return flat

def generate_subsets(pca, df_cat, ener, clustering):
    """
    Assumes df_cat is sorted in the same way as pca
    ! WARNING: assumes 6 clusters for now? rm modify colordict accordingly.

    Two input arguments: 

    pca matrix (2D)
    packing list (1D)?
    energy/dens/sg data (2D) 
    clustering labels (1D)
    """
    colordict = { 0 : "mediumturquoise", 1 : "mediumpurple" , 2: "steelblue" , 3: "hotpink" , 4: "mediumorchid", 
                  5 : "gold" , 6 : "lightpink", 7 : "darkblue", 8 :"grey", 9: "magenta", 10:"darkolivegreen", 
                  11: "lightgreen"}
    unique_clusters = df_cat.packing.nunique()
    d_clusters = {'motif':[],
                  'pca0':[],
                  'pca1' : [],
                  'color' : [],
                  'clustering_label':[], 
                  'dma_dens': [], 
                  'dma_relat_ener': [], 
                  'f_dens' : [],
                  'f_relat_ener' : []}

    df_cat = df_cat.reset_index()
    for i, r in df_cat.iterrows():
        d_clusters['motif'].append(int(r.packing))
        d_clusters['pca0'].append(pca.components_[0][i])
        d_clusters['pca1'].append(pca.components_[1][i])
        d_clusters['color'].append(colordict[int(clustering.labels_[i])])
        d_clusters['clustering_label'].append(clustering.labels_[i])
        d_clusters['dma_relat_ener'].append(ener.iloc[i]['dma_relat_ener'])
        d_clusters['dma_dens'].append(ener.iloc[i]['dma_dens'])
        d_clusters['f_dens'].append(ener.iloc[i]['f_dens'])
        d_clusters['f_relat_ener'].append(ener.iloc[i]['f_relat_ener'])
    
    d = pd.DataFrame(d_clusters)

    markers = ['o', 'v', '^', '>', 's', 'd', '*', 'X', '.'] #'>', 's','d', 'h', 'X','.']
    s0 = d[d['motif'] == 0 ]
    s1 = d[d['motif'] == 1 ]
    s2 = d[d['motif'] == 2 ]
    s3 = d[d['motif'] == 3 ]
    s4 = d[d['motif'] == 4 ]
    s5 = d[d['motif'] == 5 ]
    s6 = d[d['motif'] == 6 ]
    s7 = d[d['motif'] == 7 ]
    s8 = d[d['motif'] == 8 ]

    return s0, s1, s2, s3, s4, s5, s6, s7, s8, markers

#%%


var_dict = { 8: "nan", 0 : "none", 1 : "translational", 2 : "back-to-back", 3 : "interlocked",
            4 : "translational/back-to-back", #4 : 
            5 : "translational/interlocked",
            6 : "back-to-back/interlocked", 7 : "translational/back-to-back/interlocked"}

# inp = float(sys.argv[1]) #espilon dbscan
kernel_name = "PM_4_6"

min_sample_size = 1
n_cut = 6 #int(sys.argv[5]) #6 #11
clustering_algo = str(sys.argv[1])
n_components = int(sys.argv[2])
n_clusters = int(sys.argv[3])
linkage_criterium = str(sys.argv[4])
plotting = "False"
# print(n_cut)
# read in data:
adj = pd.read_csv(f'data/{kernel_name}.csv', header=None) # kernel
filns = pd.read_csv('data/names_list.txt', header=None, names=['ID'])
feature = pd.read_csv(f'data/Motifs601.csv')
result = pd.concat([filns,feature ], axis=1)
ener = pd.read_csv(f'data/Cyano_CrystE_all.csv')

# sort motifs in order of ID listing in Kernel:
feature = feature.set_index('ID') # make ID the index
feature = feature.reindex(filns.ID.values.tolist())#feature.loc[filns.ID.values.tolist()] # sort rows by index according to filns
feature['packing'] = feature['packing'].fillna(8) #replace(nan,8
feature['packing'].astype('int32')
feature.drop(['Unnamed: 0'], axis=1)
# add index (cyrstal_ID) to index and columns
adj = adj.set_index([filns.ID.values])
adj.columns = [filns.ID]

ener = ener.set_index('ID')
ener = ener.loc[filns.ID.values.tolist()]
labels_list = ener.sg.values.tolist()
ener['sg_label'] = ener['sg'].apply(labels_list.index)
ener['dma_relat_ener'] = ener['dma_ener'] - ener['dma_ener'].min()
ener['f_relat_ener'] = ener['f_ener'] #- ener['f_ener'].min()

# delete outlier of 100% similarity? 
# figure out how to shrink the margin between ''outlier'' similarity of graph
# with itself compared to any other crystal to better see the other similarities.
# adj = adj.replace(15057.281250, 200)
flat = see_five_highest_sim(adj)
real_max = flat[n_cut]
# print(flat[:15])
# print("Old maximum Sim: ", np.max(flat))
# print("New maximum sim: ", real_max)

#descending_x =  -np.sort(-x)
duplicates = [item for item in flat[:n_cut]]
for item in duplicates:
    adj[adj == item] = real_max #1000 #real_max #1000#real_max #np.nan # real_max #flat[-2]

flat = see_five_highest_sim(adj)
# print(flat[:15])

x = adj.values
#%%
# normalise data
norm = (x-np.min(x))/(np.max(x)-np.min(x))
G = nx.Graph(norm)
#%%
pca = PCA(n_components=n_components)
pca.fit(norm)  
#%%
X = np.transpose(pca.components_)

#%%
# Check version
if clustering_algo == "Agglomerative":
    clustering = AgglomerativeClustering(linkage=linkage_criterium, n_clusters=n_clusters, distance_threshold=None, compute_distances=True) # ward #affinity='precomputed',
elif clustering_algo == "DBSCAN":
    clustering = DBSCAN(eps=inp, min_samples=min_sample_size)
elif clustering_algo == "kmeans":
    #km = KMeans(n_clusters = 5, random_state = 1).fit(X)
    clustering = KMeans(n_clusters=n_clusters, random_state=0)

X_fitted = clustering.fit(X)

if clustering_algo == "Agglomerative":
    # print("Min cluster dist: ", min(list(X_fitted.distances_)), "max: ",max(X_fitted.distances_))
    # print("Mean Distance between clusters: ", mean(list(X_fitted.distances_)))
    print(n_clusters, min(list(X_fitted.distances_)), max(X_fitted.distances_), mean(list(X_fitted.distances_)), np.std(list(X_fitted.distances_)))

if clustering_algo == "kmeans":
    # calculate distances between cluster centres
    # distance between clusters: 
    dists = euclidean_distances(X_fitted.cluster_centers_)
    tri_dists = dists[np.triu_indices(n_clusters, 1)]
    max_dist, avg_dist, min_dist = tri_dists.max(), tri_dists.mean(), tri_dists.min()
    std = np.std(list(dists))
    print(n_clusters, min_dist, max_dist, avg_dist, std)

clusters = list(set(clustering.labels_))
n_clusters = len(set(list(clustering.labels_)))
# print('N clusters: ', n_clusters)
# #, "Distance between clusters: ", avg_dist)
if plotting == "False":
    quit()

#%%

qualitative_colors = sns.color_palette('hls', n_clusters, as_cmap=True)
qualitative_colors2 = sns.color_palette("Set3", n_clusters, as_cmap=True)
markers = ['o', 'v', '^','>', 's','d', 'h', 'X','.']

s0, s1, s2, s3, s4, s5, s6, s7, s8, m = generate_subsets(pca, feature, ener, clustering)

fig, ax = plt.subplots(figsize=(11.69,8.27)) #figure(figsize=(8.27, 11.69), dpi=100)
ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="Undefined")
ax.scatter(s1.pca0, s1.pca1,marker=m[1], c=s1.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T")
ax.scatter(s2.pca0, s2.pca1,marker=m[2], c=s2.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="B")
ax.scatter(s3.pca0, s3.pca1,marker=m[3], c=s3.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="I")
ax.scatter(s4.pca0, s4.pca1,marker=m[4], c=s4.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T,B" ) 
ax.scatter(s5.pca0, s5.pca1,marker=m[5], c=s5.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T,I")
ax.scatter(s6.pca0, s6.pca1,marker=m[6], c=s6.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="B,I" ) 
ax.scatter(s7.pca0, s7.pca1,marker=m[7], c=s7.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T,B,I")
ax.scatter(s8.pca0, s8.pca1,marker=m[8], c=s8.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="Undefined")

clusters_d = defaultdict(list)
texts= []
for idx, name in enumerate(pca.components_[1]):
    yval = name
    xval = pca.components_[0][idx]
    W99_label = adj.index[idx]  # ID label
    E_ranked = ener[ener.index == int(W99_label)]
    E_rank_label = E_ranked.iloc[0]['E_idx']
    clusters_d[clustering.labels_[idx]].append(W99_label)

    if E_rank_label <= 38: # E_rank_label energy ranking label
            texts.append(plt.text(xval,yval, E_rank_label, weight="bold", fontsize=10)) # 

    # if name < -0.10:  #and xval > 0.08:
    #     texts.append(plt.text(xval,yval, W99_label, weight="bold", fontsize=6))
    # # if name < -0.10 or name > 0.1 or xval > 0.072: #name > 0.12 or
    # if name < 0.00 and xval > 0.08:
    #     texts.append(plt.text(xval,yval, W99_label, weight="bold", fontsize=6))
    #         # print(W99_label)
    # if name > 0.3: 
    #     print("Needs to be removed before clustering: ")
    #     print(W99_label,"Idx in PCA0/PCA1: ", idx)
    #     print(x[78][90])

adjust_text(texts)
plt.xlabel("PC1", fontsize=14)
plt.ylabel("PC2", fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
        #   fancybox=True, shadow=True, ncol=4)
# plt.tight_layout()
plt.savefig(f"{kernel_name}_{n_clusters}_PCA_{linkage_criterium}_{clustering_algo}_motifs.png", dpi=300)
# plt.show()

if clustering_algo == "Agglomerative":
#importing linkage and denrogram from scipy
    #creating dendrogram based on the dummy data with single linkage criterion
    fig, ax = plt.subplots(figsize=(11.69,8.27))
    #figure(figsize=(8.27, 11.69), dpi=100)
    den = dendrogram(linkage(X, method=linkage_criterium), labels=adj.index)
    ax.set_ylabel('Distance', fontsize = 14)
    ax.set_title('Agglomerative Clustering: dendogram of 2,15-dicyano[6]helicene data', fontsize=12)
    plt.savefig(f"Dendogram_{kernel_name}_{n_clusters}_{linkage_criterium}.png",dpi=100)
    # plt.show()

#%%
# fig, ax = plt.subplots(figsize=(11.69,8.27)) #figure(figsize=(8.27, 11.69), dpi=100)
# ax.scatter(s0.dma_dens, s0.dma_relat_ener,marker=m[0], c=s0.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="Undefined")
# ax.scatter(s1.dma_dens, s1.dma_relat_ener,marker=m[1], c=s1.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T")
# ax.scatter(s2.dma_dens, s2.dma_relat_ener,marker=m[2], c=s2.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="B")
# ax.scatter(s3.dma_dens, s3.dma_relat_ener,marker=m[3], c=s3.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="I")
# ax.scatter(s4.dma_dens, s4.dma_relat_ener,marker=m[4], c=s4.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T,B" ) 
# ax.scatter(s5.dma_dens, s5.dma_relat_ener,marker=m[5], c=s5.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T,I")
# ax.scatter(s6.dma_dens, s6.dma_relat_ener,marker=m[6], c=s6.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="B,I" ) 
# ax.scatter(s7.dma_dens, s7.dma_relat_ener,marker=m[7], c=s7.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="T,B,I")
# ax.scatter(s8.dma_dens, s8.dma_relat_ener,marker=m[8], c=s8.color, s=35, alpha=0.75, cmap=qualitative_colors2, label="Undefined")

# # plt.xlabel("PC1", fontsize=14)
# # plt.ylabel("PC2", fontsize=14)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
# plt.ylim(0,25.5)
# plt.xlim(1.22, 1.32)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
#           fancybox=True, shadow=True, ncol=4)
# plt.tight_layout()
# plt.savefig(f"{kernel_name}_{n_clusters}_W99_{linkage_criterium}_{clustering_algo}_motifs601.png", dpi=300)
# # plt.show()

quit()

# DFT E-dens plot with packings and clusters: 
fig, ax = plt.subplots() #figure(figsize=(8.27, 11.69), dpi=100)
ax.scatter(s0.f_dens, s0.f_relat_ener,marker=m[0], c=s0.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="Undefined")
ax.scatter(s1.f_dens, s1.f_relat_ener,marker=m[1], c=s1.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="T")
ax.scatter(s2.f_dens, s2.f_relat_ener,marker=m[2], c=s2.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="B")
ax.scatter(s3.f_dens, s3.f_relat_ener,marker=m[3], c=s3.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="I")
ax.scatter(s4.f_dens, s4.f_relat_ener,marker=m[4], c=s4.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="T,B" ) 
ax.scatter(s5.f_dens, s5.f_relat_ener,marker=m[5], c=s5.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="T,I")
ax.scatter(s6.f_dens, s6.f_relat_ener,marker=m[6], c=s6.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="B,I" ) 
ax.scatter(s7.f_dens, s7.f_relat_ener,marker=m[7], c=s7.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="T,B,I")
ax.scatter(s8.f_dens, s8.f_relat_ener,marker=m[8], c=s8.color, s=50, alpha=0.75, cmap=qualitative_colors2, label="Undefined")

clusters_d = defaultdict(list)
texts= []
for idx, name in enumerate(pca.components_[1]):
    yval = ener.iloc[idx]['f_relat_ener'] 
    xval = ener.iloc[idx]['f_dens']
    W99_label = adj.index[idx]  # ID label
    E_ranked = ener[ener.index == int(W99_label)]
    E_rank_label = E_ranked.iloc[0]['E_idx']
    clusters_d[clustering.labels_[idx]].append(W99_label)

    if E_rank_label <= 38: # E_rank_label energy ranking label
            texts.append(plt.text(xval,yval, W99_label, weight="bold", fontsize=10)) # 

adjust_text(texts)
plt.xlabel("Density (g cm$^{-3}$)", fontsize=14)
plt.ylabel("Relative DFT-XDM energy (kJ mol$^{-1}$)", fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
# plt.ylim(0,15.5)
# plt.xlim(1.22, 1.32)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=4)
plt.tight_layout()
plt.savefig(f"{kernel_name}_{n_clusters}_DFT_{clustering_algo}_{linkage_criterium}_motifs_newK.png", dpi=300)
# plt.show()


# clusters = clustering.labels_
# clus = np.unique(clusters)
# mot = np.unique(encoded)

# stats_motif = np.zeros([6,6])

# for i in clus:
#     vals = encoded[clusters==i]
#     for j in mot:       
        
#         stats_motif[i,j-1] = len(vals[vals==j])


# # normalise by number in that motif        
# for i in clus:
#     stats_motif[:,i] = stats_motif[:,i]/np.sum(stats_motif[:,i])
      
# # normalised by number of samples in cluster
# for i in clus:
#     stats_motif[i,:] = stats_motif[i,:]/np.sum(stats_motif[i,:])
    
    
# holej = data.holeJ.values
# elecj = data.elecJ.values

# stats_holej = np.zeros([6])
# stats_elecj = np.zeros([6])


# for i in clus:
#     stats_holej[i] = np.mean(holej[clusters==i])
#     stats_elecj[i] = np.mean(elecj[clusters==i])


# # normal   

# """
# Stats is a matrix where each row describes a cluster and each column describes the number
# of a particular type of motif in that cluster (normalised to total number of that motif).

# e.g. we find that motif 1 ('herringbone') is predominantly found in cluster 0 and cluster 1. 

# We find that motif 4 (translational) is predominantly found in cluster 4

# etc etc


# """


# colors = mean_dist
# plt.scatter(X[:,0],X[:,1],c=colors)
# plt.show()

# #%% colour by min dist

# colors = min_dist
# plt.scatter(X[:,0],X[:,1],c=colors)
# plt.show()


