for i in {2..11}; do
    python motifs_clustered_copy.py $i "ward"

done

# linkages to try: ward, complete, single
# {‘ward’, ‘complete’, ‘average’, ‘single’}#
