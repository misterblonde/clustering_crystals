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


#%%

n_clusters = 6
kernel_name = "PM_4_6"
# read in data:
adj = pd.read_csv(f'data/{kernel_name}.csv', header=None) # kernel
filns = pd.read_csv('data/names_list.txt', header=None, names=['ID'])
ener = pd.read_csv(f'data/Ethynyl_CrystE_all.csv')

# add index (cyrstal_ID) to index and columns
adj = adj.set_index([filns.ID.values])
adj.columns = [filns.ID]

print(filns.head())
print(ener.head())
# crystal ener-dens file sort in order of index/col ID labels in adjacency matrix A:
ener = ener.set_index('ID')
ener = ener.loc[filns.ID.values.tolist()]
labels_list = ener.sg.values.tolist()
ener['sg_label'] = ener['sg'].apply(labels_list.index)
print(ener.head())

# delete outlier of 100% similarity? 
# figure out how to shrink the margin between ''outlier'' similarity of graph
# with itself compared to any other crystal to better see the other similarities.
# adj = adj.replace(15057.281250, 200)
x = adj.values

descending_x = -np.sort(-x)
# print(descending_x[1])
# sim_1 = []
# sim_2 = []
# for idx, col in enumerate(adj.columns.values):
# #     print(col)
# #     largest_two = adj[col].nlargest(2)
# #     sim_1.append(largest_two[0])
# #     sim_2.append(largest_two[1])
# #     print(max(sim_2))

# print("Largest similarities with itself:")
# print(min(sim_1), max(sim_1))
# print("highest similarity with something else:")
# print(min(sim_2), max(sim_2))

adj[adj > 6000] = 300
# adj.replace({15057.281250: 200}, regex=True)
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
clustering = AgglomerativeClustering(linkage='complete', n_clusters=n_clusters)
X_fitted = clustering.fit(X)
clusters = list(set(clustering.labels_))
print(n_clusters)
    #plot_clustering(X, clustering.labels_, "%s linkage" % 'complete')
    #plt.show()
# plt.figure()
# plt.scatter(X[:,0],X[:,1],c=clustering.labels_,s=20)
# plt.show()
# vmeasure.append(v_measure_score(encoded,clustering.labels_))    
# plt.show()

#%%

# right_side = ax.spines["right"]
# right_side.set_visible(False)

qualitative_colors = sns.color_palette("Set3", 8, as_cmap=True)
qualitative_colors2 = sns.color_palette("Set3", int(len(set(labels_list))), as_cmap=True)
# cmap = sns.cubehelix_palette(as_cmap=True)
# pca = PCA(n_components=2)
# pca.fit(norm)  
# n_clusters = list(set(clustering.labels_))
fig, ax = plt.subplots() # clustering.labels_
#plt.scatter(pca.components_[0],pca.components_[1], c=clustering.labels_, s=30, cmap=qualitative_colors)
# plt.scatter(pca.components_[0],pca.components_[1], c=ener.dens.values, s=30, cmap='viridis')
plt.scatter(pca.components_[0],pca.components_[1], c=ener.sg_label.values, s=30, cmap=qualitative_colors2)
# sg Pastel2_
texts = []
for idx, name in enumerate(pca.components_[1]):
    print(name)
    if name > -0.10: #or yval > 0.6:
        yval = name
        xval = pca.components_[0][idx]
        label = adj.index[idx]
        if label <= 38:
            texts.append(plt.text(xval,yval, label, weight="bold", fontsize=6)) # texts.append())

adjust_text(texts)
plt.legend(labels=labels_list)
plt.xlabel("PC1", fontsize=14)
plt.ylabel("PC2", fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
# plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='on')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.colorbar()
plt.tight_layout()
plt.savefig(f"{kernel_name}_sg_PCA.png", dpi=300)
plt.show()

quit()

#%%

clusters = clustering.labels_
clus = np.unique(clusters)
mot = np.unique(encoded)

stats_motif = np.zeros([6,6])

for i in clus:
    vals = encoded[clusters==i]
    for j in mot:       
        
        stats_motif[i,j-1] = len(vals[vals==j])


# normalise by number in that motif        
for i in clus:
    stats_motif[:,i] = stats_motif[:,i]/np.sum(stats_motif[:,i])
      
# normalised by number of samples in cluster
for i in clus:
    stats_motif[i,:] = stats_motif[i,:]/np.sum(stats_motif[i,:])
    
   
    
    
    
holej = data.holeJ.values
elecj = data.elecJ.values

stats_holej = np.zeros([6])
stats_elecj = np.zeros([6])


for i in clus:
    stats_holej[i] = np.mean(holej[clusters==i])
    stats_elecj[i] = np.mean(elecj[clusters==i])


# normal   

"""
Stats is a matrix where each row describes a cluster and each column describes the number
of a particular type of motif in that cluster (normalised to total number of that motif).

e.g. we find that motif 1 ('herringbone') is predominantly found in cluster 0 and cluster 1. 

We find that motif 4 (translational) is predominantly found in cluster 4

etc etc


"""

#%% colour by mean dist

colors = mean_dist
plt.scatter(X[:,0],X[:,1],c=colors)
plt.show()

#%% colour by min dist

colors = min_dist
plt.scatter(X[:,0],X[:,1],c=colors)
plt.show()

#%% motif type

colors = encoded
plt.scatter(X[:,0],X[:,1],c=colors,s=20)
plt.show()



