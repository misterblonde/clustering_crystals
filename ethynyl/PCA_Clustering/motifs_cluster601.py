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


def generate_subsets(pca, df_cat, clustering):
    """
    Assumes df_cat is sorted in the same way as pca
    """
    print(df_cat.packing.unique())
    unique_clusters = df_cat.packing.nunique()
    d_clusters = {'motif':[],
                  'pca0':[],
                  'pca1' : [],
                #   'pca2' : [],
                  'clustering_label':[]}

    # print(np.shape(pca.components_))
    print(np.shape(pca.components_))
    df_cat = df_cat.reset_index()
    for i, r in df_cat.iterrows():
        # print(i, r.values)
        # if i < 50:
        d_clusters['motif'].append(int(r.packing))
        d_clusters['pca0'].append(pca.components_[0][i])
        d_clusters['pca1'].append(pca.components_[1][i])
        # d_clusters['pca2'].append(pca.components_[2][i])
        d_clusters['clustering_label'].append(clustering.labels_[i])

    
    # print(d_clusters)
    d = pd.DataFrame(d_clusters)# pd.DataFrame.from_dict(d_clusters)

    qualitative_colors = sns.color_palette("Set3", 9, as_cmap=True)
    markers = ['o', 'v', '^','>','s','d', '*', 'X','.']#'>', 's','d', 'h', 'X','.']
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


# n_clusters = 6
kernel_name = "PM_4_6"
clustering_algo = "DBSCAN"
# read in data:
adj = pd.read_csv(f'data/{kernel_name}.csv', header=None) # kernel
filns = pd.read_csv('data/names_list.txt', header=None, names=['ID'])
feature = pd.read_csv(f'data/Motifs601.csv')
result = pd.concat([filns,feature ], axis=1)

# sort motifs in order of ID listing in Kernel:
feature = feature.set_index('ID') # make ID the index
feature = feature.loc[filns.ID.values.tolist()] # sort rows by index according to filns
feature['packing'] = feature['packing'].fillna(8) #replace(nan,8
feature['packing'].astype('int32')
feature.drop(['Unnamed: 0'], axis=1)
# add index (cyrstal_ID) to index and columns
adj = adj.set_index([filns.ID.values])
adj.columns = [filns.ID]

# figure out how to shrink the margin between ''outlier'' similarity of graph
# with itself compared to any other crystal to better see the other similarities.
# adj = adj.replace(15057.281250, 200)

print(adj.head())
print(feature.head())
x = adj.values

sim_1 = []
sim_2 = []
# print(adj.columns.tolist())
# print(adj.index.tolist())
# for idx, col in enumerate(adj.index.tolist()):
#     print(idx, col)
#     largest_two = adj[col,:].nlargest(2).values
#     sim_1.append(largest_two[0])
#     sim_2.append(largest_two[1])
#     print(max(sim_2))

# print("Largest similarities with itself:")
# print(min(sim_1), max(sim_1))
# print("highest similarity with something else:")
# print(min(sim_2), max(sim_2))

#adj = adj.loc[adj < 600]
# 600 = 200
print(adj.max())
adj[adj > 2000] = 600
print(adj.max())
# adj.replace({15057.281250: 200}, regex=True)
#%%
# normalise data
norm = (x-np.min(x))/(np.max(x)-np.min(x))

G = nx.Graph(norm)



#%%
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
pca.fit(norm)  


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
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans

if clustering_algo == "Agglomerative":
    clustering = AgglomerativeClustering(linkage='complete', n_clusters=n_clusters)
elif clustering_algo == "DBSCAN":
    clustering = DBSCAN(eps=1, min_samples=1)
elif clustering_algo == "kmeans":
    clustering = KMeans(n_clusters=n_clusters, random_state=0)


X_fitted = clustering.fit(X)
clusters = list(set(clustering.labels_))
n_clusters = len(set(list(clustering.labels_)))
print('N clusters: ', n_clusters)
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

qualitative_colors = sns.color_palette('hls', 9, as_cmap=True)#"Set3", 9, as_cmap=True)
markers = ['o', 'v', '^','>', 's','d', 'h', 'X','.']

s0, s1, s2, s3, s4, s5, s6, s7, s8, m = generate_subsets(pca, feature, clustering)
# pca = PCA(n_components=2)
# pca.fit(norm)  
# n_clusters = list(set(clustering.labels_))
fig, ax = plt.subplots()
# plt.scatter(pca.components_[0],pca.components_[1], marker=markers[feature.packing.values], c=clustering.labels_, s=30, cmap=qualitative_colors, label=clustering.labels_)
ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
ax.scatter(s1.pca0, s1.pca1,marker=m[1], c=s1.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
ax.scatter(s2.pca0, s2.pca1,marker=m[2], c=s2.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )

ax.scatter(s3.pca0, s3.pca1,marker=m[3], c=s3.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
ax.scatter(s4.pca0, s4.pca1,marker=m[4], c=s4.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ ) #ax

ax.scatter(s5.pca0, s5.pca1,marker=m[5], c=s5.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
ax.scatter(s6.pca0, s6.pca1,marker=m[6], c=s6.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ ) #ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
ax.scatter(s7.pca0, s7.pca1,marker=m[7], c=s7.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
ax.scatter(s8.pca0, s8.pca1,marker=m[8], c=s8.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )

# annotate scatter points with crystal nomer
# texts = []
# for idx, name in enumerate(pca.components_[1]):
#     print(name)
#     if name > -0.10: #or yval > 0.6:
#         yval = name
#         xval = pca.components_[0][idx]
#         label = adj.index[idx]
#         texts.append(plt.text(xval,yval, label, weight="bold", fontsize=6)) # texts.append())

# adjust_text(texts)
plt.xlabel("PC1", fontsize=14)
plt.ylabel("PC2", fontsize=14)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.colorbar(plt.cm.ScalarMappable(cmap=qualitative_colors), ax=ax)
plt.tight_layout()
plt.savefig(f"{kernel_name}_{n_clusters}_PCA_{clustering_algo}_motifs601.png", dpi=300)
plt.show()

# fig, ax = plt.subplots()
# pca1 = pca2
# # plt.scatter(pca.components_[0],pca.components_[1], marker=markers[feature.packing.values], c=clustering.labels_, s=30, cmap=qualitative_colors, label=clustering.labels_)
# ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s1.pca0, s1.pca1,marker=m[1], c=s1.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s2.pca0, s2.pca1,marker=m[2], c=s2.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s5.pca0, s5.pca1,marker=m[3], c=s5.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s6.pca0, s6.pca1,marker=m[4], c=s6.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ ) #ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s7.pca0, s7.pca1,marker=m[5], c=s7.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )
# ax.scatter(s8.pca0, s8.pca1,marker=m[6], c=s8.clustering_label, s=60, cmap=qualitative_colors, label=clustering.labels_ )

# # annotate scatter points with crystal nomer
# # texts = []
# # for idx, name in enumerate(pca.components_[1]):
# #     print(name)
# #     if name > -0.10: #or yval > 0.6:
# #         yval = name
# #         xval = pca.components_[0][idx]
# #         label = adj.index[idx]
# #         texts.append(plt.text(xval,yval, label, weight="bold", fontsize=6)) # texts.append())

# # adjust_text(texts)
# plt.xlabel("PC2", fontsize=14)
# plt.ylabel("PC3", fontsize=14)
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.colorbar(plt.cm.ScalarMappable(cmap=qualitative_colors), ax=ax)
# plt.tight_layout()
# plt.savefig(f"{kernel_name}_{n_clusters}_PCA23_{clustering_algo}_motifs601.png", dpi=300)
# plt.show()


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


colors = mean_dist
plt.scatter(X[:,0],X[:,1],c=colors)
plt.show()

#%% colour by min dist

colors = min_dist
plt.scatter(X[:,0],X[:,1],c=colors)
plt.show()


