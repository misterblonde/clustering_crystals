for i in {2..12}; do
    python motifs_clustered.py $i "ward"

done

# linkages to try: ward, complete, single
# {‘ward’, ‘complete’, ‘average’, ‘single’}#
