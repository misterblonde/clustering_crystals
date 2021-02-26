# PM 4 6 kernel only!
# GRID SEARCH

# test KMeans clluster sizes and PCA dimensions:
for i in {2..11}; do
    python motifs_clustered_copy.py "kmeans" 2 $i "ward" >>Kmeans_PCA2.txt
    python motifs_clustered_copy.py "kmeans" 3 $i "ward" >>Kmeans_PCA3.txt
    python motifs_clustered_copy.py "kmeans" 6 $i "ward" >>Kmeans_PCA6.txt

done

# linkages to try: ward, complete, single
# {‘ward’, ‘complete’, ‘average’, ‘single’}#
# agglomerative clusterings:
for i in {2..11}; do
    python motifs_clustered_copy.py "Agglomerative" 2 $i "ward" >>Agglo_PCA2_ward.txt
    python motifs_clustered_copy.py "Agglomerative" 3 $i "ward" >>Agglo_PCA3_ward.txt
    python motifs_clustered_copy.py "Agglomerative" 6 $i "complete" >>Agglo_PCA6_ward.txt

    python motifs_clustered_copy.py "Agglomerative" 2 $i "complete" >>Agglo_PCA2_complete.txt
    python motifs_clustered_copy.py "Agglomerative" 3 $i "complete" >>Agglo_PCA3_complete.txt
    python motifs_clustered_copy.py "Agglomerative" 6 $i "complete" >>Agglo_PCA6_complete.txt

    python motifs_clustered_copy.py "Agglomerative" 2 $i "single" >>Agglo_PCA2_single.txt
    python motifs_clustered_copy.py "Agglomerative" 3 $i "single" >>Agglo_PCA3_single.txt
    python motifs_clustered_copy.py "Agglomerative" 6 $i "single" >>Agglo_PCA6_single.txt
done
