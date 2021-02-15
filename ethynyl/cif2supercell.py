from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.core.structure import Structure
import numpy as np
import sys
import pywindow as pw
import pandas as pd
# from pymatgen.analysis.graphs import StructureGraph
# from pymatgen.analysis.local_env import MinimumDistanceNN
from scipy.spatial import distance
import subprocess as sp
# import pytest
"""

PYMATGEN ALTERNATIVE TO BUILDING THE SUPERCELL FROM CIF FILE

requirements:
* ! only works in the py36 env
- pymatgen
- pywindow (python 3.6, git clone from repos)

ideas to improve this:
- get the neighbouring dimers according to axis directions (?)

INPUT:
- cif file
- n_atoms (single molecule)

OUTPUT:
- dimer xyz (nonreal)

"""


def find_supercell_centre(basen):
    """
    takes the supercell (XYZ) and finds the centre-of-mass (COM)
    of the supercell.

    make sure supercell size is large enough to get every dimer
    created, without any missing atoms!

    Input:
        - basen --input argument: filn of supercell.
    returns:
        - com: centre-of-mass of supercell
        - coms_list: list of all centre-of-mass of all mols in supercell
        - all_mols: list of xyz of all molecules in supercell.

    This function is the only part of the script requiring pywindow.
    """
    molsys = pw.MolecularSystem.load_file(f"{basen}_supercell.xyz")
    # read supercell
    molsys.make_modular()  # finds separate molecules in supercell

    center_center = []
    # list of all the COMs of the mols in supercell (locations.txt)
    all_mols = []
    for molecule in molsys.molecules:
        mol = molsys.molecules[molecule]
        info = mol.calculate_centre_of_mass()
        # previously: full_analysis() get all molecule information in
        # pymatgen access info.keys()
        if mol.no_of_atoms == n_atoms:  # if molecule is not broken append to
            center_center.append(info.tolist())  # com to locations list
            all_mols.append(mol)  # molecule to mols list

    # find centroid of all COMs from supercell
    center_center = list(center_center)
    # print("Center of supercell found.")
    coms = len(center_center)
    print(f"There are {coms} molecules with the correct #atoms in supercell.")

    assert coms != 0, "Error: Incorrect n_atoms per mol entered"

    # find centre of locations (coms)
    x = [row[0] for row in center_center]
    y = [row[1] for row in center_center]
    z = [row[2] for row in center_center]
    centroid = (sum(x)/coms, sum(y)/coms, sum(z)/coms)

    return centroid, center_center, all_mols


def closest_node(node, nodes):
    """
    Function finds the molecule closest to the centre of mass of the supercell.
    (closest distance by COM-COM distance)

    returns
        - closest molecule to COM of supercell
        - idx of closest molecule to COM of supercell.

    """
    closest_index = distance.cdist([node], nodes).argmin()

    return nodes[closest_index], closest_index


def duplicate_dimers_finder(dist_only):
    """
    Function checks whether current distance is same as previous distance

    Returns
        - idx of all duplicate positions

    * TODO: there are easier way to do this in pandas, like keepfirst,
            * but didn't work.
    """
    dist_only = dist_only.round(6)  # avoid differences due to numerical noise
    idx_duplicates = []
    last_item = None
    for i, dist in dist_only.iteritems():
        # print("i: dist ", i, dist)
        if dist == last_item:
            idx_duplicates.append(i)
            # print("duplicate found")
        last_item = dist  # needs to be inside if statement!

        # print(idx_duplicates)

    return idx_duplicates


def find_closest_packed_dimers(closest, closest_idx, coms_list):
    """

    We only require the closest dimers, not dimers which are further apart
    so in reality they would have another molecule sitting in between.

    Hence we sort all mols by their distance from our most central molecule
    Then, we only select the 5 closest unique dimers, ignoring dimers that
    are duplicates ie. have the same COM-COM distance.

    returns
        - 5 unique nearest neighbour dimers

    """
    # find all mols closest to the one closest to the centroid of the supercell
    df = pd.DataFrame(np.array(coms_list), columns=['x', 'y', 'z'])
    df["dist"] = np.sqrt(((df.x-closest[0])**2) +
                        ((df.y-closest[1])**2)+((df.z-closest[2])**2))

    # idx of df indicates molecules to be used/at closest distance
    sorted_df = df.sort_values(by=["dist"])

    dist_only = pd.Series(sorted_df['dist'], dtype="float64")
    idx_duplicates = duplicate_dimers_finder(dist_only)

    # let's not remove duplicates for graph clustering of crystals,
    # thus keep lowest 14 molecules:
    clean = sorted_df
    # remove duplicate dimer motifs with the same COM-COM distance:
    # clean = sorted_df.drop(idx_duplicates)


    # save unique positions to file (optional)
    sorted_df.to_csv("coms_distance_sorted.csv")
    return clean


def generate_sm_xyz(central_top, all_mols):
    # Generates single molecule real_xyz files
    for idx in central_top:
        all_mols[idx].dump_molecule(
            "./mol_{0}.pdb".format(idx),
            include_coms=False,
            override=True)
        sp.call(f"babel -ipdb mol_{idx}.pdb -oxyz mol_{idx}.xyz", shell=True)
    sp.call("rm mol_*.pdb", shell=True)


def generate_dimer_xyz(idx_closest, central_top):
    """

    Function generates real dimer xyz files and opens them.

    """
    # remove header from central molecule file
    sp.call(
        f"tail -n +3 mol_{idx_closest}.xyz > nonreal_{idx_closest}.xyz",
        shell=True)


    # make neighshell file from monomers

    # makes dimer files from monomers
    for idx, name in enumerate(central_top[1:]):  # start from 2nd mol in list
        sp.call(f"tail -n +3 mol_{name}.xyz > nonreal_{name}.xyz", shell=True)


    # make neighshell minRep (minimum representation )
    sp.call(
        f"cat header_neighshell nonreal_*.xyz > \
            real_minRep_{basen}.xyz", shell=True)


    sp.call(
        f"cat nonreal_*.xyz > \
            minRep_{basen}.xyz", shell=True)
    #
    #     # open dimers in vesta
    #     sp.call(f"open real_dimer_{idx_closest}_{name}.xyz", shell=True)


#
# OUTPUT settings:
out_mol = 14
# number of unique nearest neighbour dimers you want to
# extract from crystal

# Input handling
filn = sys.argv[1]
basen = filn.split(".")[0]
n_atoms = int(sys.argv[2])  # make sure to input correct thing here
parser = CifParser(filn)
struc = parser.get_structures()[0]


s = Structure.from_file(filn)
scaling_matrix = [[4, 0, 0], [0, 4, 0], [0, 0, 4]] # [[4, 0, 0], [0, 4, 0], [0, 0, 4]]  # create 4x4x4 matrix
s.make_supercell(scaling_matrix, to_unit_cell=False)
print(np.shape(s))  # n
xyzrep = XYZ(s)
xyzrep.write_file(f"{basen}_supercell.xyz")  # write supercell to file

# find centre-of-mass of supercell and return all COMs and mols xyz.
centroid, coms_list, all_mols = find_supercell_centre(basen)

# find molecule closest to centroid of supercell
closest, idx_closest = closest_node(centroid, coms_list)

# get list of all the unique closest dimers in ascending order
clean = find_closest_packed_dimers(closest, idx_closest, coms_list)

# make dimer xyz files for (25 or) 15 closest molecules
central_top = clean[:out_mol+1].index.values.tolist()

# makes single molecule xyz files (rm pdb files created during this)
generate_sm_xyz(central_top, all_mols)

# make dimer xyz
generate_dimer_xyz(idx_closest, central_top)
# use the closest 5-6 mols as dimers
