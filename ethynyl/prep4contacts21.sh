


# * SCRIPT: This script prepares the dimers to be in the correct format
# * for using the ShortContactsCode. 
#  Input:   dir_name , must contain dir_name/real_*xyz 
#           (nonreal xyz files)
#
#
#  Output: dir inside ShortContacts folder, containing:
#           1)  .mol files of dimers 
#           2)  names.txt base_name list of dimers in folder
#

d=$1


cd $d

# loop through files in folder
for file in *.xyz; do
    base_name=$(basename $file ".xyz")
    
    if [ "$base_name" == "sm" ];then
        :
    else
        echo $base_name >> names.txt
        
        # number of atoms in dimer
        natoms=$(wc -l $file | awk '{print $1}')
        
        
        babel -ixyz $file -omol "$base_name".mol
        
        #initial_folder=$PWD
        #cd "../../3_short_contacts_code/"
        
        #if [ ! -d $d ]; then
        #    mkdir $d
        #fi
        
        #cp "$inital_folder"/"$base_name".mol $d/
        
    fi
    
done

#mv names.txt $d/
cd ../

# copy the newly generated .mol files to 3_short_contacts_code for analysis
#cp -r $d/ ../../../3_short_contacts_code/"$d"/

