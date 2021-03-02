k="PM_4_6"
for i in {3..6}; do
    python motifs_clustered.py "kmeans" 2 $i "ward" 13 $k            #>>Kmeans_"$k"_PCA2.txt
    python motifs_clustered.py "Agglomerative" 2 $i "ward" 13 $k     #>>Kmeans_"$k"_PCA2_ward.txt
    python motifs_clustered.py "Agglomerative" 2 $i "complete" 13 $k #>>Kmeans_"$k"_PCA2_complete.txt
    # python motifs_clustered_copy.py "kmeans" 3 $i "ward" 13 $k >>Kmeans_"$k"_PCA3.txt
    # python motifs_clustered_copy.py "kmeans" 6 $i "ward" 13 $k >>Kmeans_"$k"_PCA6.txt

done

k="PM_5_7"
for i in {3..6}; do
    python motifs_clustered.py "kmeans" 2 $i "ward" 13 $k            #>>Kmeans_"$k"_PCA2.txt
    python motifs_clustered.py "Agglomerative" 2 $i "ward" 13 $k     #>>Kmeans_"$k"_PCA2_ward.txt
    python motifs_clustered.py "Agglomerative" 2 $i "complete" 13 $k #>>Kmeans_"$k"_PCA2_complete.txt
    # python motifs_clustered_copy.py "kmeans" 3 $i "ward" 13 $k >>Kmeans_"$k"_PCA3.txt
    # python motifs_clustered_copy.py "kmeans" 6 $i "ward" 13 $k >>Kmeans_"$k"_PCA6.txt

done
