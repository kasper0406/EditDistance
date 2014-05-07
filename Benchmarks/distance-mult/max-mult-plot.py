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
plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,2))

mult_rand = pd.read_csv("max-multiply-random.dat", delim_whitespace=True).head(30)
mult_short = pd.read_csv("max-multiply-short.dat", delim_whitespace=True).head(30)
mult_long = pd.read_csv("max-multiply-long.dat", delim_whitespace=True).head(30)
slow_mult = pd.read_csv("slow-max-multiply.dat", delim_whitespace=True).head(8)

ax1.plot(mult_rand['N'], mult_rand['median'] / mult_rand['N'], label="Random shuffle")
ax1.plot(mult_short['N'], mult_short['median'] / mult_short['N'], label="Short candidate list")
ax1.plot(mult_long['N'], mult_long['median'] / mult_long['N'], label="Long candidate list")
ax1.plot(slow_mult['N'], slow_mult['median'] / slow_mult['N'], label="Naive running time")

#ax2 = ax1.twinx()
#ax2.set_yscale('log', basey=2)
#ax2.plot(mult_long['N'], mult_long['instructions'] / mult_long['N'], label="Instructions", ls='dashed', color='c')

ax1.legend()
plt.savefig('max-dist-mult.pdf')
