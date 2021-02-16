
# runs entire workflow on one crystal
n=$1 # crystal no. 
n_atoms=46

source activate py27

cd opt_res/ 

bash convert.sh $n

cp opt_"$n".cif ../ 

cd ../

source deactivate py27 

source activate py36 

bash run-dimer-Ebind.sh opt_"$n".cif $n_atoms 

#bash prep4contacts21.sh opt_"$n"

conda deactivate 

