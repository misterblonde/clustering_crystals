#!/bin/bash


inp_folder=$1
n_atoms=$2
out_folder=$3

mkdir $out_folder

python NeighShell_from_Supercells.py $inp_folder $n_atoms $out_folder


exit 0
