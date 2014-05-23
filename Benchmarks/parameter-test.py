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

input_data = pd.read_csv("parameter-test-random.dat", delim_whitespace=True)

xstep, xlimit = 0.1, 2
fstep, flimit = 0.5, 5

data = [[0 for x in xrange(int(flimit / fstep) - 1)] for x in xrange(int(xlimit / xstep) - 1)]
for i in input_data.index:
  data[int(input_data.xfactor[i] * (1 / xstep) - 1)][int(input_data.ffactor[i] * (1 / fstep) - 1)] = input_data.median_time[i]
  # data[int(input_data.xfactor[i] - 1)][int(input_data.ffactor[i] - 1)] = input_data.time[i]
data = np.array(data)

fig, ax = plt.subplots(1)
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

cbar = fig.colorbar(heatmap)
cbar.set_label("Running time [s]")

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(0, flimit / fstep, 1) + 0.5, minor=False)
ax.set_yticks(np.arange(0, xlimit / xstep, 1) + 0.5, minor=False)

xlabels = np.around([(x + 1) * xstep for x in xrange(int(xlimit / xstep) - 1)], 2)
flabels = np.around([(x + 1) * fstep for x in xrange(int(flimit / fstep) - 1)], 2)

ax.set_xticklabels(flabels, minor=False)
ax.set_yticklabels(xlabels, minor=False)

ax.set_ylim([0, int(xlimit / xstep) - 1])
ax.set_xlim([0, int(flimit / fstep) - 1])

ax.set_title("Block size parameters")
ax.set_xlabel("ffactor")
ax.set_ylabel("xfactor")

plt.savefig('parameter-test-random.pdf')




input_data = pd.read_csv("parameter-test-fib.dat", delim_whitespace=True)

xstep, xlimit = 0.3, 6.3
fstep, flimit = 0.3, 6.3

data = [[0 for x in xrange(int(flimit / fstep) - 1)] for x in xrange(int(xlimit / xstep) - 1)]
for i in input_data.index:
  data[int(input_data.xfactor[i] * (1 / xstep) - 1)][int(input_data.ffactor[i] * (1 / fstep) - 1)] = input_data.median_time[i]
  # data[int(input_data.xfactor[i] - 1)][int(input_data.ffactor[i] - 1)] = input_data.time[i]
data = np.array(data)

fig, ax = plt.subplots(1)
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

cbar = fig.colorbar(heatmap)
cbar.set_label("Running time [s]")

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(0, flimit / fstep, 1) + 0.5, minor=False)
ax.set_yticks(np.arange(0, xlimit / xstep, 1) + 0.5, minor=False)

xlabels = np.around([(x + 1) * xstep for x in xrange(int(xlimit / xstep) - 1)], 2)
flabels = np.around([(x + 1) * fstep for x in xrange(int(flimit / fstep) - 1)], 2)

ax.set_xticklabels(flabels, minor=False)
ax.set_yticklabels(xlabels, minor=False)

ax.set_ylim([0, int(xlimit / xstep) - 1])
ax.set_xlim([0, int(flimit / fstep) - 1])

ax.set_title("Block size parameters")
ax.set_xlabel("ffactor")
ax.set_ylabel("xfactor")

plt.savefig('parameter-test-fib.pdf')
