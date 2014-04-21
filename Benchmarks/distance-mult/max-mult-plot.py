import sys
sys.path.insert(0, "/Library/Python/2.7/site-packages")

import matplotlib as mpl
# mpl.use("pgf")
# pgf_with_pdflatex = {
#   "pgf.texsystem": "pdflatex",
#   "pgf.preamble": [
#     r"\usepackage[utf8x]{inputenc}",
#     r"\usepackage[T1]{fontenc}",
#     r"\usepackage{cmbright}"
#   ]
# }
# mpl.rcParams.update(pgf_with_pdflatex)

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

plt.xkcd()

# Min distance product base case plot
plt.figure()
plt.title("Varying base case size")
ax1 = plt.subplot(111)
plt.xscale('log', basex=2)
#plt.yscale('log', basey=2)
plt.ylabel('Normalized running time')
plt.xlabel('N')

mult = pd.read_csv("max-multiply.dat", delim_whitespace=True)
slow_mult = pd.read_csv("slow-max-multiply.dat", delim_whitespace=True).head(7)

ax1.plot(mult['N'], mult['median'] / mult['N'], label="Algorithm running time")
ax1.plot(slow_mult['N'], slow_mult['median'] / slow_mult['N'], label="Naive running time")

# ax2 = ax1.twinx()
# ax2.set_yscale('log', basey=2)
# ax2.plot(mult['N'], mult['L3_miss'] / mult['N'], label="Instructions", ls='dashed', color='k')

ax1.legend()
plt.savefig('max-dist-mult.pdf')
