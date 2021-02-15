import community
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import argparse
# Replace this with your networkx graph loading depending on your format !

if __name__ == '__main__':

    # Pass Computed Kernel
    parser = argparse.ArgumentParser(description='A test file for all kernels')

    parser.add_argument('-i','--SimilarityKernel', \
                                    help='Similarity Matrix from Grakel Kernel')
    args = parser.parse_args()

    df = pd.read_csv(args.SimilarityKernel, index_col=0)

    print(df.head())
    results = pd.read_csv('results.csv')
    print(results.head())

    G = nx.from_numpy_matrix(df.values)



    # Compute the best partition (Louvain Method)
    partition = community.best_partition(G)

    # Plot
    size = float(len(set(partition.values())))
    pos = nx.spring_layout(G)
    count = 0.
    for com in set(partition.values()) :
        count = count + 1.
        list_nodes = [nodes for nodes in partition.keys()
                                    if partition[nodes] == com]

        print(str(count/size))
        #G.nodes()
        # results = results.set_index('name')
        # results = results.reindex(G.nodes())
        #
        # results['motif']=pd.Categorical(results['motif'])
        # results['motif'].cat.codes

        nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 15,
    #    node_color=results['motif'].cat.codes, cmap=plt.cm.Set1,) #20
                                    node_color =str(count / size), cmap=plt.cm.viridis) #node_color

    print(count, list_nodes)
    #nx.draw_networkx_labels(G, pos)   # postition in matrix
        #colors = {'D':'red', 'E':'blue', 'F':'green', 'G':'black'}
    #nx.draw_networkx_edges(G, pos, edge_color ='magenta', edge_size=1, alpha=0.5)
    plt.show()



# color coding
# -120 - 170 range of energies if used as color label
#     node_color = []
#     for node in g.nodes(data=True):
# # QUESTION:
#
#     # if the node has the attribute group1
#         if 'group1' in node[1]['group']:
#             node_color.append('blue')
#
#         # if the node has the attribute group1
#         elif 'group2' in node[1]['group']:
#             node_color.append('red')
#
#         # if the node has the attribute group1
#         elif 'group3' in node[1]['group']:
#             node_color.append('green')
#
#         # if the node has the attribute group1
#         elif 'group4' in node[1]['group']:
#             node_color.append('yellow')
#
#         # if the node has the attribute group1
#         elif 'group5' in node[1]['group']:
#             node_color.append('orange')
