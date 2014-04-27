import sys
sys.path.insert(0, "/Library/Python/2.7/site-packages")

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

plt.xkcd()
plt.xlabel('N')
plt.ylabel('Normalized time')

df = pd.read_csv("lcs_blowup_total.dat", delim_whitespace=True)

# normalized_times = df['median'] / (df['N'] * ((2 * df['N']) / (df['A_prod'] + df['B_prod'])) * np.sqrt(np.log2((2 * df['N']) / (df['A_prod'] + df['B_prod']))))
normalized_times = 10E3 * df['median'] / (df['N'] * (df['A_prod'] + df['B_prod']) * np.sqrt(np.log2(df['N'] / (df['A_prod'] + df['B_prod']))))

# plt.ylim([0,0.1])
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)
plt.plot(df['N'], normalized_times, label="Runningtime")

#for (name, df_type) in df.groupby('type'):
#  if name == "total":
#    continue;
#
#  aggregated = df_type.groupby('N').agg([ np.mean, np.max, np.min ]).reset_index()
#  normalized_times = aggregated['time']['mean'] / aggregated['N']

#  plt.plot(aggregated['N'], normalized_times, label=name)

plt.legend()
plt.show()
