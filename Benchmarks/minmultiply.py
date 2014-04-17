import sys
sys.path.insert(0, "/Library/Python/2.7/site-packages")

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

plt.xkcd()
plt.xlabel('n')
plt.ylabel('Normalized time')

df = pd.read_csv("minmultiply.dat", sep='\t')

for (name, df_type) in df.groupby('type'):
  if name == "total":
    continue;

  aggregated = df_type.groupby('n').agg([ np.mean, np.max, np.min ]).reset_index()
  normalized_times = aggregated['time']['mean'] / aggregated['n']

  plt.plot(aggregated['n'], normalized_times, label=name)

plt.legend()
plt.show()
