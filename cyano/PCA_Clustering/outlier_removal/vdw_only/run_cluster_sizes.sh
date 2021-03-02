# PM 4 6 kernel only!
# GRID SEARCH
k=$1 # kernel

# test KMeans clluster sizes and PCA dimensions:
for i in {2..11}; do
    python motifs_clustered_copy.py "kmeans" 2 $i "ward" 13 $k >>Kmeans_"$k"_PCA2.txt
    python motifs_clustered_copy.py "kmeans" 3 $i "ward" 13 $k >>Kmeans_"$k"_PCA3.txt
    python motifs_clustered_copy.py "kmeans" 6 $i "ward" 13 $k >>Kmeans_"$k"_PCA6.txt

done

# linkages to try: ward, complete, single
# {‘ward’, ‘complete’, ‘average’, ‘single’}#
# agglomerative clusterings:
for i in {2..11}; do
    python motifs_clustered_copy.py "Agglomerative" 2 $i "ward" 13 $k >>Agglo_"$k"_PCA2_ward.txt
    python motifs_clustered_copy.py "Agglomerative" 3 $i "ward" 13 $k >>Agglo_"$k"_PCA3_ward.txt
    python motifs_clustered_copy.py "Agglomerative" 6 $i "ward" 13 $k >>Agglo_"$k"_PCA6_ward.txt

    python motifs_clustered_copy.py "Agglomerative" 2 $i "complete" 13 $k >>Agglo_"$k"_PCA2_complete.txt
    python motifs_clustered_copy.py "Agglomerative" 3 $i "complete" 13 $k >>Agglo_"$k"_PCA3_complete.txt
    python motifs_clustered_copy.py "Agglomerative" 6 $i "complete" 13 $k >>Agglo_"$k"_PCA6_complete.txt

    python motifs_clustered_copy.py "Agglomerative" 2 $i "single" 13 $k >>Agglo_"$k"_PCA2_single.txt
    python motifs_clustered_copy.py "Agglomerative" 3 $i "single" 13 $k >>Agglo_"$k"_PCA3_single.txt
    python motifs_clustered_copy.py "Agglomerative" 6 $i "single" 13 $k >>Agglo_"$k"PCA6_single.txt
done

mkdir -p GridSearch_"$k"

mv Kmeans_"$k"*.txt GridSearch_"$k"
mv Agglo_"$k"*.txt GridSearch_"$k"
