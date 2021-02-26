for file in *.txt; do
    echo $file

    python plot_grid_search.py $file
done
