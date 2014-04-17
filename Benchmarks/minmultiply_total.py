import sys
sys.path.insert(0, "/Library/Python/2.7/site-packages")

from math import log
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

plt.xkcd()
plt.xlabel('n')
plt.ylabel('Time')

# Simple algorithm
simple = pd.read_csv("minmultiply_simple.dat", sep='\t')
aggregated = simple.groupby('n').agg([ np.median, np.max, np.min ]).reset_index()
# normalized_times = aggregated['time']['median'] / (aggregated['n'] * log(aggregated['n'], 2))
plt.plot(aggregated['n'], aggregated['time']['median'], label="Simple")

# Divide and Conquer algorithm
d_and_q = pd.read_csv("minmultiply.dat", sep='\t')
total = d_and_q.groupby('type').get_group('total')
aggregated = total.groupby('n').agg([ np.median, np.max, np.min ]).reset_index()
# normalized_times = aggregated['time']['median'] / (aggregated['n'] * log(aggregated['n'], 2))
plt.plot(aggregated['n'], aggregated['time']['median'], label="Divide & Conquer (bc = 1)")

# Divide and Conquer (bc = 20) algorithm
d_and_q = pd.read_csv("minmultiply_dandq.dat", sep='\t')
total = d_and_q.groupby('type').get_group('total')
aggregated = total.groupby('n').agg([ np.median, np.max, np.min ]).reset_index()
# normalized_times = aggregated['time']['median'] / (aggregated['n'] * log(aggregated['n'], 2))
plt.plot(aggregated['n'], aggregated['time']['median'], label="Divide & Conquer (bc = 20)")

plt.legend()
plt.show()
