#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

ANALYSIS OF GRAPH SIMILARITY KERNELS

1. PCA plot
2. Clustering
@author: JAS

"""

def see_five_highest_sim(adj):
    x = adj.values
    flat = x.flatten()
    flat = list(set(flat))
    flat.sort(reverse=True)
    print(flat[:5])
    return flat


def check_av_n_mol_common(clusters_d):
    """
    For each crystal in respective cluster, check how many mols it has in common
    with other crystals within the cluster. 
    """
    df = pd.read_csv('data/n_mols_comp.csv', header=None, names=['id1', 'sg1','id2', 'sg2', 'n_matched', 'rmsd'])
    print("We are checking the clusters now:")
    df.id1.astype("int32")
    df.id2.astype("int32")
    for k,v in clusters_d.items():
        print(k,v)
        for i in v:
            for j in v:
                if i != j:
                    print(i, j) # works fine up till here 
                    mols_in_cluster = df[(df['id1'] == int(i)) & (df['id2'] == int(j))] 
                    # if mols_in_cluster.empty() == False: 
                    print(mols_in_cluster)


    
#%%
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns 
from adjustText import adjust_text


#%%

n_clusters = 6
kernel_name = "PM_4_6"
# read in data:
adj = pd.read_csv(f'data/{kernel_name}.csv', header=None) # kernel
filns = pd.read_csv('data/names_list.txt', header=None, names=['ID'])
ener = pd.read_csv(f'data/Cyano_CrystE_all.csv')

# add index (cyrstal_ID) to index and columns
adj = adj.set_index([filns.ID.values])
adj.columns = [filns.ID]

# crystal ener-dens file sort in order of index/col ID labels in adjacency matrix A:
ener = ener.set_index('ID')
ener = ener.loc[filns.ID.values.tolist()]
labels_list = ener.sg.values.tolist()
ener['sg_label'] = ener['sg'].apply(labels_list.index)
ener['dma_relat_ener'] = ener['dma_ener'] - ener['dma_ener'].min()

x = adj.values
# delete outlier of 100% similarity? 
# figure out how to shrink the margin between ''outlier'' similarity of graph
# with itself compared to any other crystal to better see the other similarities.
# adj = adj.replace(15057.281250, 200)
print("Make sure you adapt the number of duplicate/highest similarity outliers \
      to rescale for clustering.")
flat = see_five_highest_sim(adj)
duplicate1 = flat[0]
duplicate2 = flat[1]
real_max = flat[4]
sim_with_itself = flat[2]
#descending_x =  -np.sort(-x)
duplicates = [flat[0], flat[1], flat[2], flat[3]]
for item in duplicates:
    adj[adj == item] = real_max #np.nan # real_max #flat[-2]


see_five_highest_sim(adj)
# %%

#%%
# normalise data
norm = (x-np.min(x))/(np.max(x)-np.min(x))
G = nx.Graph(norm)

#%%

# data = pd.read_csv('results.csv')
# min_dist = data.min_dist.values
# mean_dist = data.mean_dist.values


# motif = pd.get_dummies(data.motif)

# for i in range(0,6):
#     motif.iloc[:,i] = motif.iloc[:,i]*(i+1)

# encoded = motif.sum(axis=1).values

#%%
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
pca.fit(norm)  
# plt.scatter(pca.components_[0],pca.components_[1])
#%%

X = np.transpose(pca.components_)

#%%
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import v_measure_score
# vmeasure = []
# for i in range(2,20):
#     clustering = AgglomerativeClustering(linkage='complete', n_clusters=i)
#     clustering.fit(norm)
#     #plot_clustering(X, clustering.labels_, "%s linkage" % 'complete')
#     #plt.show()
#     plt.figure()
#     plt.scatter(norm[:,0],norm[:,1],c=clustering.labels_,s=20)
 
    # print(v_measure_score)
    # vmeasure.append(v_measure_score(encoded,clustering.labels_))
    # plt.show()
#%%
clustering = AgglomerativeClustering(linkage='ward', n_clusters=n_clusters)
X_fitted = clustering.fit(X)
clusters = list(set(clustering.labels_))
print("N clusters", n_clusters)
# print(clustering.labels_)
#%%
qualitative_colors = sns.color_palette("Set3", 8, as_cmap=True)
qualitative_colors2 = sns.color_palette("Set3", int(len(set(labels_list))), as_cmap=True)
# cmap = sns.cubehelix_palette(as_cmap=True)
# pca = PCA(n_components=2)
# pca.fit(norm)  
# n_clusters = list(set(clustering.labels_))
fig, ax = plt.subplots() # clustering.labels_
#plt.scatter(pca.components_[0],pca.components_[1], c=clustering.labels_, s=30, cmap=qualitative_colors)
# plt.scatter(pca.components_[0],pca.components_[1], c=ener.dens.values, s=30, cmap='viridis')
plt.scatter(pca.components_[0],pca.components_[1], c=clustering.labels_, s=30, cmap=qualitative_colors2)
# sg Pastel2_
texts = []

from collections import defaultdict
clusters_d = defaultdict(list)

for idx, name in enumerate(pca.components_[1]):
    # print(name)
    # if name < -0.40 or name > 0.40: #or yval > 0.6:
        # print("PCA1 val: ", name, "index label: ", adj.index[idx])
    yval = name
    xval = pca.components_[0][idx]
    W99_label = adj.index[idx]  # ID label
    E_ranked = ener[ener.index == int(W99_label)]
    E_rank_label = E_ranked.iloc[0]['E_idx']
    print(W99_label, E_rank_label, clustering.labels_[idx])

    clusters_d[clustering.labels_[idx]].append(W99_label)
    # print(E_ranked)

    # if E_rank_label <= 38: # E_rank_label energy ranking label
    #         texts.append(plt.text(xval,yval, label, weight="bold", fontsize=6)) # 
    if name > 0.12 or name < -0.05 or xval > 0.12:
            texts.append(plt.text(xval,yval, W99_label, weight="bold", fontsize=6))

adjust_text(texts)
# plt.legend(labels=labels_list)
plt.xlabel("PC1", fontsize=14)
plt.ylabel("PC2", fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
# plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.colorbar()
plt.tight_layout()
plt.savefig(f"{kernel_name}_{n_clusters}_sg_PCA_newKernel.png", dpi=300)
plt.show()


# check_av_n_mol_common(clusters_d)

# quit()

# #%%

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

# #%% colour by mean dist

# colors = mean_dist
# plt.scatter(X[:,0],X[:,1],c=colors)
# plt.show()

# #%% colour by min dist

# colors = min_dist
# plt.scatter(X[:,0],X[:,1],c=colors)
# plt.show()

# #%% motif type

# colors = encoded
# plt.scatter(X[:,0],X[:,1],c=colors,s=20)
# plt.show()



