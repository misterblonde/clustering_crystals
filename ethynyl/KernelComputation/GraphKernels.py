import numpy as np
import pandas as pd
import os, re, sys
import grakel
import time


# INPUT DATA
current_path = "../" #os.getcwd()
#path = "/Users/julia/Documents/GitHub/Crystals2Graphs/csv_NEIGHSHELL/"
inp_folder = "csv" #"csv_NEIGHSHELL_test"

# KERNEL SELECTION
req_kernel = "pm" #str(sys.argv[1]) #"nh"


# LABELS GENERATOR
n_atoms = 44 # number of atoms in single molecule
shell_size = 15  # number of neighbouring moleucles in shell
c_end = 28



def KernelSimilarity(new_graphs, node_labels):
    """

    Combines Node Labels and adds them onto graphs which are read in as .csv files
    Builds a similarity matrix which is then saved as a .csv file

    """

    graphs = []
    for ref_index, ref_file in enumerate(new_graphs):

        # Compute Reference Kernel all other Kernels will be compared across the row:
        # Generate graph that all others are compared to in this row:
        df = pd.read_csv(f"{current_path}/{inp_folder}/{ref_file}.csv", header=None)
        G = df.values.tolist()
        # join networkX graph data and fake labels created previously:
        graphs.append([G, node_labels])


    start = time.time()
    out = kernel.fit_transform(graphs)
    end = time.time()
    print("Time needed for Reference Kernel computation (s): ", end - start)

    # Write Matrix to file
    pd.DataFrame(out, index=new_graphs, columns=new_graphs).to_csv(f"{req_kernel}_{inp_folder}_SimilarityKernel.csv", header=True, index=True)



# node labels assuming cyanomolecule
def label_maker(n_atoms, c_end, n ):
    # GENERATE NODE LABELS for grakel graph representation
    node_labels = {}
    atom_no =0
    for n_mols in range(0, shell_size):
        #print(n_mols)
        for i in range(0,n_atoms):
            #for element in node_labels.values():
            if i < c_end:
                node_labels[atom_no]= 'C'
            elif i< c_end+2:
                node_labels[atom_no] = 'N'
            else:
                node_labels[atom_no] = 'H'
            #print(i)
            atom_no += 1
            #print(atom_no)

    return node_labels


node_labels = label_maker(n_atoms, c_end, n=2)


# kernels sp kernel for now
if req_kernel == "sp": # or "shortest_path":
    kernel = grakel.GraphKernel(kernel={"name": "shortest_path"}, labels=None)
elif req_kernel == "pm": # or "pyramid_match":
    kernel = grakel.GraphKernel(kernel={"name": "pyramid_match"}, normalize=True, verbose=True) # with_labels=True, L=2
elif req_kernel == "rw": # or "random_walk":
    kernel = grakel.GraphKernel(kernel={"name": "random_walk"}, normalize=True)
elif req_kernel == "nh": # or "neighborhood_hash":
    kernel = grakel.GraphKernel(kernel={"name": "neighborhood_hash"}, normalize=True)
else:
    print("Error: The kernel you're requesting doesn't exist.")
    print("Possible kernels are: shortest_path (sp), pyramid_match (pm), random_walk (rw), neighborhood_hash (nh)")
    exit()


print(f"The selected kernel parameters are {kernel.get_params()}")
##############################################################################
# GRAKEL KERNELS TRIED
# fast
# random walk - values range 0 - (-1)
# pyramid pyramid_match - all vals close to 1
# ***neighborhood_hash - good range of no. 0-1

# slow
# shortest Path
# lovasz_theta

# not managed
# random walk labeled
# weisfeiler_lehman - not a base kernel

# modify class attributes, doesn't change the outcome, but speed
#pm_kernel.L = 2
#pm_kernel.L = 12
#kernel=pm_kernel
##############################################################################

# rigid list of all files in dir
comp_crystals = os.listdir(f"{current_path}/{inp_folder}")

comp_crystals
# strip .csv extension from files in directory
new = [x.strip(".csv") for x in comp_crystals]

# output DataFrame
df_ =pd.DataFrame(index=new, columns=new)
#############################################
train_graphs = []
for ref_index, ref_file in enumerate(new):
    # Generate graph that all others are compared to in this row:
    # df = pd.read_csv(f"{current_path}/{inp_folder}/{ref_file}.csv", header=None)
    # G = df.values.tolist()
    # graph = [[G, node_labels]]
    # train_graphs.append(graph)
    
    G, node_labels = read_gpickle(f"{current_path}/{inp_folder}/{ref_file}.pkl")

    grakel_graph = grakel.graph_from_networkx([G], node_labels_tag=node_labels)
#train_graphs

#print(np.shape(train_graphs))

# train set fit_transform
ref_kernel = kernel.fit_transform(train_graphs)
#ref_kernel
print(ref_kernel)

print(np.shape(ref_kernel))


df_.loc[ref_file, ref_file]= ref_kernel.flat[0]
print(f"{ref_file}\t{ref_file}\t{ref_kernel.flat[0]}")
end = time. time()
print("Time needed for Reference Kernel computation (s): ", end - start)
# reference kernel if normalised always 1.0

#comp_crystals
for index, name in enumerate(new):
    #print(name)
# generate all other similarity scores with respect to reference one of row:
    print("Current file is: ", name)
    start2 = time. time()
    graph_file = pd.read_csv(f"{current_path}/{inp_folder}/{name}.csv", header=None)#dtype=int,
    graph_list= graph_file.values.tolist()
    grakel_graph = [[graph_list, node_labels]]
    #print("Shape: ", np.shape(grakel_graph), "\n")
    #print("grakel graph\n", grakel_graph)
    #print("col 2 grakel graph\n", grakel_graph[:,1])

    # test set _transform
    similarity = kernel.transform(grakel_graph)
    print(f"{ref_file}\t{name}\t{similarity.flat[0]}")
    end2 = time. time()
    print("Time needed for Kernel comparison (s): ", end2 - start2)

##############################################################################

    # write similarity to output similarity matrix
        df_.loc[ref_file,name] = similarity.flat[0]

#
df_.to_csv(f"{req_kernel}_{inp_folder}_SimilarityKernel.csv", header=True, index=True)
#####################################################################
