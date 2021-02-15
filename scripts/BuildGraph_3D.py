import numpy as np
import sys
import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import subprocess as sp

"""       CLUSTERING STEP II

        BUILD GRAPH FROM AN XYZ 3x3 SUPERCELL

        Tried for res3_supercell - gives one connected graph (VERSION 1)
        However if res1_supercell - graph is disconnected (VERSION 1)

        VERSION 2: res3 and res1 are connected



        GRAPH DISTANCE THRESHOLD:

        VERSION 1: only distances ≤ sum of the VdW radii will be considered

        VERSION 2.2: Distances ≤ sum of the VdW radii + 2 Angstrom will be considered

"""


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

def text_to_int(text,char2idx):
    return np.array([char2idx[c] for c in text])

class simple_binary_graph:

    def __init__(self,name, labels, coords):
        self.name = name
        self.labels = labels
        self.coords = coords
        self.dist_mat = distance_matrix(coords, coords)
        self.graph_dim = len(labels)


    def vdw_distance_graph(self):
        # turns distance matrix into a simple 0 1 binary graph according to
        # sum of VdW distance + 5 A.
        # get vdw radii for each atom to compute sum of vdw radii
        #sp.call("mkdir csv", shell=True)
        labels_vdw= []
        for i in range(0,len(self.labels)):
            vdw = atom_vdw_radii.get(str(self.labels[i]))
            labels_vdw.append(vdw)

        # turning distance matrix into an unweighted graph - simplest version
        # decide what kind of distance matarix you want:
        counter=0
        for i in range(0,int(np.shape(self.labels)[0])):
            for j in range(0,int(np.shape(self.labels)[0])):
                sum_vdw_radii = labels_vdw[i]+labels_vdw[j]
                #print("{} vs. {}".format(sum_vdw_radii, dist_mat[i, j]))
                threshold = float(sum_vdw_radii) + 2
                if float(self.dist_mat[i, j]) <= threshold: #+3
                    #print("{} <= {}".format( self.dist_mat[i, j], sum_vdw_radii))
                #if float(dist_mat[i, j]) <= sum_vdw_radii:
                    self.dist_mat[i,j]= 1
                #else:
                    #dist_mat[i,j]=(dist_mat[i,j]-sum_vdw_radii)/(threshold-sum_vdw_radii)

                    counter =counter + 1

                else:
                    self.dist_mat[i,j]= 0


        G = nx.from_numpy_matrix(self.dist_mat)
        final_mat = nx.to_scipy_sparse_matrix(G)#.toarray() # (84,84)

        labels_arr = np.zeros((self.graph_dim,self.graph_dim))
        for idx_i, i in enumerate(self.labels):
            for idx_j, j in enumerate(self.labels):

                new_str = i+j
                #print(new_str)
                if new_str == "CC":
                    bond_type = 1
                elif (new_str == "CH") or (new_str == "HC"):
                    bond_type = .5
                elif new_str == "HH":
                    bond_type = 0
                else:
                    print("Unknown bond type.")

                labels_arr[idx_i][idx_j] = bond_type


        A = final_mat.toarray() # graph Aij
        B = labels_arr # atom labels

        A = A.reshape(1,self.graph_dim,-1)
        B = B.reshape(1,self.graph_dim,self.graph_dim)

        final3D = np.dstack([A, B])
        final3D = np.reshape(final3D, (2,self.graph_dim,self.graph_dim))

        final2D = pd.DataFrame(A.toarray())
        final2D.to_csv(f"./csv/{self.name}.csv", index=None, header=None)

        with open(f"./csv/{self.name}.pkl",'wb') as f:
            pickle.dump(final3D, f)
        

        return self.dist_mat


def atom_type():
    return
################ SUPERCELL SPECIFIC VARIABLES #############################
inp_folder = str(sys.argv[1])
# input folder name
# input folder should contain numbers of the pairs eg dimer_0_1.xyz

import glob, re
def main():
    #names = []
    for name in glob.glob("./{}/*xyz".format(inp_folder)):
        #names.append(name)
        #print(name)
        try:
            basename = re.findall("\d+_\d+_\d+", name)[0]
            #print(basename)
            labels = np.genfromtxt(name,usecols=0,dtype=str, skip_header=0)
            coords = np.genfromtxt(name, skip_header=0)[:,1:]

            graph = simple_binary_graph(basename, labels, coords)
            binary_dist_mat = graph.vdw_distance_graph()

            # Graph from distance matrix
            G = nx.from_numpy_matrix(binary_dist_mat)#, labels=graph_labels)
            #better_mat = nx.to_scipy_sparse_matrix(G)
            #print(type(better_mat))
            # check if there are no loose atoms anywhere:
            if nx.is_connected(G)== False:
                print("Warning! {} is disconnected".format(basename))

            # number of connected components within graph
            elif nx.number_connected_components(G) > 1:
                print("Two disconnected components.")

            else:
                pass
        except:
            print(f"{name} cannot be converted into a graph with depth2.")


if __name__=="__main__":
    main()
