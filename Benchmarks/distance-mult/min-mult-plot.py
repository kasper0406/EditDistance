import sys
sys.path.insert(0, "/Library/Python/2.7/site-packages")

import matplotlib as mpl
mpl.use("pgf")
pgf_with_pdflatex = {
  "pgf.texsystem": "pdflatex",
  "pgf.preamble": [
    r"\usepackage[utf8x]{inputenc}",
    r"\usepackage[T1]{fontenc}"
  ]
}
mpl.rcParams.update(pgf_with_pdflatex)

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

# Min distance product base case plot
plt.figure()
plt.title("Varying base case size")
ax = plt.subplot(111)
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)
plt.ylabel('Normalized running time')
plt.xlabel('N')

mult = pd.read_csv("min-multiply-bc.dat", delim_whitespace=True)

colors = [ 'k', 'b', 'g', 'r', 'c', 'm', 'y' ]
sizes = [ 1, 17, 14, 20, 10, 24, 50 ]
i = 0
for (bcs, df) in mult.groupby('bcs'):
  if bcs not in sizes:
    continue

  plt.plot(df['N'], df['median'] / (df['N'] * np.log2(df['N'])), label=str(bcs), color=colors[i])

  i = i + 1

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.85])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 4)
plt.savefig('min-dist-mult-bc.pdf')

# bcs = 50 examination
plt.figure()
plt.title("Examination of base case size 50")
ax1 = plt.subplot(111)
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)
plt.ylabel('Normalized running time')
plt.xlabel('N')

ax2 = ax1.twinx()
ax2.set_ylabel('Approximate base case instance size')

df = mult[(mult.bcs == 50) & (mult.N > 50)]

plot1 = ax1.plot(df['N'], df['median'] / (df['N'] * np.log2(df['N'])), label="Normalized running time", color='b')
# plt.plot(df['N'], df['instructions'] / (df['N'] * np.log2(df['N'])), label="Instructions")
plot2 = ax2.plot(df['N'], df['N'] / (np.power(2, np.ceil(np.log2(df['N'] / 50)))), label="Approximate instance base case size", color='r')

for tl in ax1.get_yticklabels():
  tl.set_color('b')
for tl in ax2.get_yticklabels():
  tl.set_color('r')

# Add labels
plots = plot1 + plot2
labels = [ l.get_label() for l in plots ]

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])
box = ax2.get_position()
ax2.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])

ax1.legend(plots, labels, loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True)

plt.savefig('min-dist-mult-bc50.pdf')