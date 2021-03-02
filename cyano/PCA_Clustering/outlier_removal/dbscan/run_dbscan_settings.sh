# PM 4 6 kernel only!
# GRID SEARCH

#  python motifs_clustered_copy.py "DBSCAN" 2 2 0.5 5
# test KMeans clluster sizes and PCA dimensions:

# python motifs_clustered_copy.py "DBSCAN" 2 0.5 5 6 >>DBSCAN_PCA2_default.txt
# python motifs_clustered_copy.py "DBSCAN" 3 0.5 5 6 >>DBSCAN_PCA3_default.txt
# python motifs_clustered_copy.py "DBSCAN" 6 0.5 5 6 >>DBSCAN_PCA6_default.txt

# epsilon
for i in {1..7}; do

    for j in {1..7}; do
        python motifs_clustered_copy.py "DBSCAN" 2 $i $j 6 >>DBSCAN_PCA2_"$i"_"$j".txt
        python motifs_clustered_copy.py "DBSCAN" 3 $i $j 6 >>DBSCAN_PCA3_"$i"_"$j".txt
        python motifs_clustered_copy.py "DBSCAN" 6 $i $j 6 >>DBSCAN_PCA6_"$i"_"$j".txt
    done
done

for i in {0..7}; do
    python motifs_clustered_copy.py "DBSCAN" 6 "1."$i 1 6 >>DBSCAN_PCA6_epsilon_"$i".txt
done
