import grakel 
import pickle 
import networkx as nx 
from grakel import Graph
from grakel.kernels import NeighborhoodSubgraphPairwiseDistance # requires edge labels 
from grakel.kernels import NeighborhoodHash, RandomWalk, GraphHopper, GraphletSampling
import time
import os
import numpy as np
"""
This is a trial case on how to compute Similarity Kernels using Grakel 
"""

def read_in_grakel(working_dir, folder, filn):
    # read graph from pickle
    G = nx.read_gpickle(f"{working_dir}/{folder}/{filn}.pkl")
    # compute and add node_labels on betweenness centrality
    node_labels = nx.betweenness_centrality(G, normalized=True)
    nx.set_node_attributes(G, node_labels, 'label')
    # turn nx graph into grakel graph
    return grakel.graph_from_networkx([G], node_labels_tag='label')


# input parameters
current_path = ".."
inp_folder = "csv"
ref_file = "126"
req_kernel = "sp"

# read graph from pickle
G = nx.read_gpickle(f"{current_path}/{inp_folder}/{ref_file}.pkl")
# compute and add node_labels on betweenness centrality
node_labels = nx.betweenness_centrality(G, normalized=True)
nx.set_node_attributes(G, node_labels, 'label')
# turn nx graph into grakel graph
grakel_graph = grakel.graph_from_networkx([G], node_labels_tag='label')

start = time.time()
# compute kernel for graph (incl. labels)
gk = GraphletSampling()#, k=20) #normalize=True)
# print(np.shape(grakel_graph))
# fit or fit_transform?
ref_kernel = gk.fit_transform(grakel_graph)

end = time. time()
print("Time needed for Reference Kernel computation (s): ", end - start)
# reference kernel if normalised always 1.0

############# OTHER CRYSTALS #########################################
# rigid list of all files in dir
comp_crystals = os.listdir(f"{current_path}/{inp_folder}")

# strip .csv extension from files in directory
new = [x.strip(".pkl") for x in comp_crystals]

mylist = []
graphs = []
#comp_crystals
print("Generating graph objects.")
# test set _transform
print("Starting Kernel computation.")
start2 = time. time()
for index, name in enumerate(new):
        grakel_graph = read_in_grakel("..", "csv", name)
        graphs.append(grakel_graph)
        mylist.append(name)
        similarity = gk.transform(grakel_graph)
        print(f"{ref_file}\t{name}\t{similarity.flat[0]}")

end2 = time. time()
print("Time needed for Kernel comparison (s): ", end2 - start2)

##############################################################################
#
similarity.to_csv(f"{req_kernel}_{inp_folder}_SimilarityKernel.csv", header=True, index=True)
#####################################################################
with open('SimilarityKernelIndex.txt', 'w') as f:
    for item in my_list:
        f.write("%s\n" % item)