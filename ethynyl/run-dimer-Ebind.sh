#"""
#
#Generates supercell and unique dimers from the inputs:
#
#
#   Inputs required: 
#  1. CIF file
#  2. # atoms per mol
#
#   OUTPUT:
#   - xyz coords of unique dimers
#"""

#inputs:
cifname=$1 # cif file (with extension)
natoms=$2 # n_atoms per mol
ndimer=$(($natoms*15))

# get name without extension
base_name=$(basename $cifname ".cif")

# if directory doesn't exist:
if [ ! -d $base_name ]; then

mkdir "$base_name"

fi

cd $base_name

# create real_xyz dimer header
# for use of make_real.sh
cat > header_neighshell << EOF
$ndimer

EOF

cp ../"$cifname" .
#cp ../make_real.sh .

# generate dimers
python ../cif2supercell.py $cifname $natoms


#bash make_real.sh


#while read row; do

# echo $row
#
# dist=
#
# done
#cp mol_0.xyz sm.xyz
#cat nonreal_0.xyz nonreal_1.xyz nonreal_2.xyz nonreal_3.xyz nonreal_4.xyz > NN4.xyz

#rm mol*.pdb
rm header_*
rm *_supercell.xyz
rm mol*.xyz
#rm dimer*.xyz
rm nonreal*.xyz
cd ../
