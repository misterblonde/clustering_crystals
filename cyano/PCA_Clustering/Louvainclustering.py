
import community
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from itertools import chain
import seaborn as sns
import re
import matplotlib.colors as colors
from matplotlib.cm import RdBu_r, Blues
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, to_rgba
from adjustText import adjust_text


def see_five_highest_sim(adj):
    x = adj.values
    flat = x.flatten()
    flat = list(set(flat))
    flat.sort(reverse=True)
    print(flat[:5])
    return flat

def generate_subsets(whole_louvain, packings):
    """
    Assumes df_cat is sorted in the same way as pca/clustering
    """
    print(packings.packing.unique())
    unique_packings = packings.packing.nunique()

    d_clusters = {'motif':[],
                  'x':[],
                  'y' : [],
                  'clustering_label':[]}

    # # print(np.shape(pca.components_))
    # print(np.shape(pca.components_))
    df_cat = packings.reset_index()


    for i, r in df_cat.iterrows():
        # print(i, r.values)
        # if i < 50:
        d_clusters['motif'].append(int(r.packing))
        d_clusters['x'].append(whole_louvain.iloc[i]['x'])
        d_clusters['y'].append(whole_louvain.iloc[i]['y'])
        # d_clusters['pca2'].append(pca.components_[2][i])
        d_clusters['clustering_label'].append(whole_louvain.iloc[i]['class'])

    
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

def plot_scatter_clusters(s0, s1, s2, s3, s4, s5, s6, s7, s8, m, E_rank, whole, n_clusters, filn):
    # SCATTERPLOT

    fig, ax = plt.subplots()
    

    qualitative_colors = sns.color_palette("hls", n_clusters, as_cmap=True)

    # plt.scatter(pca.components_[0],pca.components_[1], marker=markers[feature.packing.values], c=clustering.labels_, s=30, cmap=qualitative_colors, label=clustering.labels_)
    ax.scatter(s0.x, s0.y,marker=m[0], c=s0.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)
    ax.scatter(s1.x, s1.y,marker=m[1], c=s1.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)
    ax.scatter(s2.x, s2.y,marker=m[2], c=s2.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)

    ax.scatter(s3.x, s3.y,marker=m[3], c=s3.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)
    ax.scatter(s4.x, s4.y,marker=m[4], c=s4.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70) #ax

    ax.scatter(s5.x, s5.y,marker=m[5], c=s5.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)
    ax.scatter(s6.x, s6.y,marker=m[6], c=s6.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70) #ax.scatter(s0.pca0, s0.pca1,marker=m[0], c=s0.clustering_label, s=30, cmap=qualitative_colors, label=clustering.labels_ )
    ax.scatter(s7.x, s7.y,marker=m[7], c=s7.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)
    ax.scatter(s8.x, s8.y,marker=m[8], c=s8.clustering_label, s=30, cmap=qualitative_colors, alpha=0.70)

    # annotate scatter points with crystal nomer
    texts = []
    
    for i, name in enumerate(E_rank.E_rank.values.tolist()):
        # print(name)
        if int(name) < 39: #or yval > 0.6:
            yval = whole.iloc[i]['y']#[0]
            xval = whole.iloc[i]['x']#[0] #pca.components_[0][idx]
            label = name #+1]\
            W99_label = E_rank.index.values[i]
            texts.append(plt.text(xval,yval, W99_label, weight="bold", fontsize=6)) # texts.append())

    adjust_text(texts)
    plt.xlabel("TSNE1", fontsize=14)
    plt.ylabel("TSNE2", fontsize=14)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # plt.colorbar(plt.cm.ScalarMappable(cmap=qualitative_colors), ax=ax)
    plt.tight_layout()
    plt.savefig(f"{outn}_{n_clusters}_{kernel_inp}.png", dpi=300)
    plt.show()


kernel_inp="PM_4_6"
outn = "Louvain"
n_cut = 13 #12

#df = pd.read_csv(f"data/{kernel_inp}.csv", header=0)
adj = pd.read_csv(f'data/{kernel_inp}.csv', header=None) # kerne
print(f"{np.shape(adj)[0]} training examples in Similarity file")

filns = pd.read_csv('data/names_list.txt', header=None, names=['ID']) # need order to cp order for feature df
feature = pd.read_csv(f'data/Motifs601.csv')
e_rank = pd.read_csv(f'data/Cyano_CrystE_all.csv')
# result = pd.concat([filns,feature], axis=1)

# sort motifs in order of ID listing in Kernel:
feature = feature.set_index('ID') # make ID the index
feature = feature.loc[filns.ID.values.tolist()] # sort rows by index according to filns
feature['packing'] = feature['packing'].fillna(8) #replace(nan,8
feature['packing'].astype('int32')
feature.drop(['Unnamed: 0'], axis=1)
print(feature.head())
# # sort erank in order of filns: 
e_rank['E_rank'] = e_rank.index
e_rank = e_rank.set_index('ID')
e_rank_id = e_rank.loc[filns.ID.values.tolist()] 
print(e_rank_id)

# add index (cyrstal_ID) to index and columns
adj = adj.set_index([filns.ID.values])
adj.columns = [filns.ID]

flat = see_five_highest_sim(adj)
real_max = flat[n_cut]
print(flat[:15])
print("Old maximum Sim: ", np.max(flat))
print("New maximum sim: ", real_max)

#descending_x =  -np.sort(-x)
duplicates = [item for item in flat[:n_cut]]
for item in duplicates:
    adj[adj == item] = real_max #1000 #real_max #1000#real_max #np.nan # real_max #flat[-2]

print(flat[:5])
x = adj.values
# figure out how to shrink the margin between ''outlier'' similarity of graph
# with itself compared to any other crystal to better see the other similarities.
# adj = adj.replace(15057.281250, 200)

# print(adj.head())
# # print(feature.head())
# adj[adj > 7000] = 0
# print(adj.max())
#x = adj.values
norm = (x-np.min(x))/(np.max(x)-np.min(x))
x = norm


#%%


# ________________ START LOUVAIN ________________________
# turn similiarty matrix into graph
G = nx.from_numpy_matrix(x)

# Louvain Method to compute best partition
#first compute the best partition
partition = community.best_partition(G)

# get clustesr
unique_clusters =set(partition.values())
print("Unique Clusters: ", unique_clusters)

# put all row_idx and class labels into a df
part = pd.DataFrame({'row':list(partition),
                     'class':list(partition.values())} )

# embedd cluster differences as two axes/cols of df
X_embedded = TSNE(n_components=2).fit_transform(x)
tsne = pd.DataFrame(data=X_embedded)
tsne = tsne.rename(columns={0: "x", 1: "y"}) #, 2:"z"})
#tsne = pd.DataFrame({'x':X_embedded[:,0], 'y':X_embedded[:,1], 'y':X_embedded[:,2]})
# print(tsne.head())
# combine both dfs
whole = pd.concat([tsne, part], axis=1)

# print(whole.head())
# n unique clusters found
n_clusters = len(unique_clusters)
# labels = np.arange(0, n_clusters, 1)

""""____________________  START PLOTTING YOUR CLUSTERS ______________"""

print("Unique clusters: ", n_clusters)


s0, s1, s2, s3, s4, s5, s6, s7, s8, m = generate_subsets(whole, feature)


print("Print features: ")

print(feature.head())
# for idx, r in s0.iterrows():
#     print(r)
    # print(r.motif[0], r.clustering_label[0])
plot_scatter_clusters(s0, s1, s2, s3, s4, s5, s6, s7, s8, m, e_rank_id, whole, n_clusters, outn)