import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import sys

filn=sys.argv[1]
basen = filn.split(".")[0]
d = pd.read_csv(filn, header=None, sep=" ", names=['n_clusters', 'min_d', 'max_d', 'av_d', 'std_d'])

print(d.head())

plt.title(f"{basen}")
plt.xlim(2,11,1)
plt.plot(d.n_clusters.values, d.max_d.values, label="maximum distance")
plt.plot(d.n_clusters.values, d.av_d.values, label="average distance")
plt.plot(d.n_clusters.values, d.std_d.values, label="std(distances)")
plt.legend()
plt.savefig(f"{basen}.png", dpi=200)
plt.show()