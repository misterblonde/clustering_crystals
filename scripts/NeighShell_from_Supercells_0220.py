#!/usr/bin/env python
import numpy as np
import pandas as pd
import re, os, sys
from pathlib import Path
import operator

################################################################################
#
#           Reduce 3x3x3 supercell to a Neighbouring Shell
#
#                       NeighShell_all.py version 0.1
#
#           Converts 3x3x3 supercell into a central molecule + 14 neighbouring
#           (=15 molecules) minimum representation of the supercell.
#
#
# Neighbouring shell is created according to 14 closes centroid centroid distance
# contacts eg. find the centre of mass of each molecule and calculate the distance
#between them, then select the lowest 14-15 molecules nearby?
###############################################################################
# INPUT VARIABLES
# out_folder = "NEIGHSHELL"    # out folder that should contain all NeighShell files
# inp_folder = "SUPERCELLS_INPUT"  # folder containing input XYZ supercells - shouldn't contain any numbers
# n_atoms = 44                # int number of atoms in single molecule
# cutoff = 15.0               # Centroid-Centroid distance cutoff in angstroms
# n_neighbours = 14
def neighbours_within_cutoff(mol, inp_folder, atoms_per_mol, cutoff):
    # Extract raw data from xyz file (i.e. discard atom symbols)
    labels = np.genfromtxt(f"{inp_folder}/{mol}",usecols=0,dtype=str)
    print(f"{int(len(labels)/int(atoms_per_mol))} molecules in supercell.")

    # init
    counter = 0
    distance_dict= {}

    for j in range(1,int(len(labels)/atoms_per_mol)):

        index_start = int(atoms_per_mol)*j
        index_end = int(atoms_per_mol)*(j+1)

        coordfile1 = np.genfromtxt(f"{inp_folder}/{mol}")[0:atoms_per_mol,1:]
        coordfile2 = np.genfromtxt(f"{inp_folder}/{mol}")[index_start:index_end,1:]

        # initialise features:
        centre1=[0,0,0]
        centre2=[0,0,0]

        MASS=0

        # Find centre of mass
        for i in range(0,atoms_per_mol):

            if labels[i]=='C':

                mass=12

            elif labels[i]=='N':

                mass=14

            elif labels[i]=='S':

                mass=32

            elif labels[i]=='P':

                mass=31

            else:

                mass=1

            MASS+=mass

            centre1+=(coordfile1[i,:]*mass)
            centre2+=(coordfile2[i,:]*mass)

        centre1/=MASS
        centre2/=MASS

        # Find distance between molecules
        distance=np.linalg.norm(centre1-centre2)

        if distance<cutoff:
            # append moleucles of interests and their distances to dict:
            distance_dict[j]=distance
            counter += 1

        return distance_dict, counter

# Input arguments ___________________________________________________________
inp_folder = str(sys.argv[1])    # folder containing input XYZ supercells
n_atoms = int(sys.argv[2])       # int number of atoms in single molecule
out_folder = str(sys.argv[3])    # outfolder which contains all NeighShell files

# Constants___________________________________________________________________
# could be passed as arguments if needed
cutoff = 15.0               # Centroid-Centroid distance cutoff in angstroms
n_neighbours = 14



for mol in os.listdir(f"{inp_folder}/"):
    # input filn
    res_name = re.findall('res([0-9.]*[0-9]+)', f"{inp_folder}/{mol}")[-1]
    # out filn
    end_file = Path(f"{out_folder}/{res_name}_minRep.xyz")


    # skip neighbouring shells that have already been created
    if end_file.is_file():
        pass
    else:
        distance_dict, counter = neighbours_within_cutoff(mol, inp_folder, n_atoms, cutoff)

        # final results for given crystal:
        print(f"Molecules within {cutoff} A cutoff: {counter}")
        #print(sorted(distance_dict.items()))
        ########################################################################

        # sort molecules within range in ascending order
        sorted_x = sorted(distance_dict.items(), key=operator.itemgetter(1))

        # select the 14 lowest distance neighbouring molecules
        final_select = []
        for mol_idx in range(0, n_neighbours):
            final_select.append(sorted_x[mol_idx])

        #########################################################################
        # I/O Write Output to file
        # idx of all molecules neighbouring central :
        mols_idx_keep = [x[0] for x in final_select]


        data = pd.DataFrame([])
        # write central molecule to file:
        init_file = pd.read_csv(f"{inp_folder}/{mol}", sep=" ", header=None)
        mol0 = init_file.iloc[0:n_atoms, :]
        data = data.append(mol0, ignore_index=True)
        #mol0.to_csv("NeighShell_{}.xyz".format(0), sep=" ", header=None, index=False)

        # write all xyz of neighbouring molecules to files:
        for z in mols_idx_keep:
            index_start = int(n_atoms)*z
            index_end = int(n_atoms)*(z+1)
            # make sure to include Atom Labels here:
            file_mol = init_file.iloc[index_start:index_end,:]
            data = data.append(file_mol, ignore_index=True)
            # write file
            #file_mol.to_csv("NeighShell_{}.xyz".format(z), sep=" ", header=None, index=False)


        data.to_csv(end_file, sep=" ", float_format='%.6f', header=None, index=False)
