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

# Fibonacci plots
plt.figure()
plt.title("Relative compression of Fibonacci strings")
ax = plt.subplot(111)
plt.xscale('log', basex=2)
plt.ylabel('Compression overhead factor')
plt.xlabel('Length')

fib = pd.read_csv("compression_fib.dat", delim_whitespace=True)
relativeProductions = (fib['compression'] * fib['N']) / (fib['num'] + 2)
normalizedRelativeProductions = relativeProductions / np.log2(fib['N'])

plt.plot(fib['N'], relativeProductions, label="Production overhead")
plt.plot(fib['N'], normalizedRelativeProductions, label="Normalized production overhead")

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.85])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True)
plt.savefig('fib.pdf')


# Random string plots
plt.figure()
plt.title("Compression quality")
plt.xscale('log', basex=2)
plt.ylabel('Compression ratio')
plt.xlabel('Length')

random = pd.read_csv("compression_random.dat", delim_whitespace=True)
# normalizedCompression = random['compression'] * (np.log2(random['N'] / 4) / np.log2(4))

plt.plot(random['N'], random['compression'], label="Random string")
# plt.plot(random['N'], normalizedCompression, label="Normalized compression")

plt.legend()
plt.savefig('random.pdf')


# DNA Plots
plt.figure()
plt.xscale('log', basex=2)
plt.ylabel('Compression ratio')
plt.xlabel('Prefix length')

genome1 = pd.read_csv("compression_fasta_genome1.dat", delim_whitespace=True)
genome2 = pd.read_csv("compression_fasta_genome2.dat", delim_whitespace=True)

plt.plot(genome1['N'], genome1['compression'], label="Genome 1")
plt.plot(genome2['N'], genome2['compression'], label="Genome 2")
plt.plot(random['N'], random['compression'], label="Random")

plt.legend()
plt.savefig('dna.pdf')


# Human genome plots
plt.figure()
plt.title("Compressibility of the human genome to a SLP")
plt.xscale('log', basex=2)
plt.ylabel('Compression ratio')
plt.xlabel('Prefix length')

hg_repetitive = pd.read_csv("compression_fasta_hg_repetitive.dat", delim_whitespace=True)
hg_nonrepetitive = pd.read_csv("compression_fasta_hg_nonrepetitive.dat", delim_whitespace=True)
hg_combined = pd.read_csv("compression_fasta_hg_combined.dat", delim_whitespace=True)

plt.plot(hg_repetitive['N'], hg_repetitive['compression'], label="HG Repetitive")
plt.plot(hg_nonrepetitive['N'], hg_nonrepetitive['compression'], label="HG Nonrepetitive")
plt.plot(hg_combined['N'], hg_combined['compression'], label="HG Combined")
plt.plot(random['N'], random['compression'], label="Random string")

plt.legend()
plt.savefig('hg.pdf')


# Running time plots
plt.figure()
ax1 = plt.subplot(111)
plt.title("Runing time for compression")
plt.xscale('log', basex=2)
ax1.set_ylabel('Normalized running time')
plt.xlabel('Length')
plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,2))

clipped_random = random.tail(14)
clipped_genome1 = genome1.tail(12)
clipped_genome2 = genome2.tail(12)
clipped_hg_repetitive = hg_repetitive.tail(14)
clipped_hg_nonrepetitive = hg_nonrepetitive.tail(14)
clipped_hg_combined = hg_combined.tail(14)
clipped_fib = fib.tail(15)

ax1.plot(clipped_fib['N'], clipped_fib['median'] / (clipped_fib['N'] * np.log2(clipped_fib['N'])), label="Fibonacci")
ax1.plot(clipped_hg_repetitive['N'], clipped_hg_repetitive['median'] / (clipped_hg_repetitive['N'] * np.log2(clipped_hg_repetitive['N'])), label="HG Repetitive")
ax1.plot(clipped_hg_nonrepetitive['N'], clipped_hg_nonrepetitive['median'] / (clipped_hg_nonrepetitive['N'] * np.log2(clipped_hg_nonrepetitive['N'])), label="HG Nonrepetitive")
ax1.plot(clipped_hg_combined['N'], clipped_hg_combined['median'] / (clipped_hg_combined['N'] * np.log2(clipped_hg_combined['N'])), label="HG Combined")
ax1.plot(clipped_random['N'], clipped_random['median'] / (clipped_random['N'] * np.log2(clipped_random['N'])), label="Random", color='k')

ax1.legend()

#ax2 = ax1.twinx()
#ax2.set_ylabel('Normalized cache misses')
# ax2.plot(clipped_fib['N'], clipped_fib['L3_miss'] / (clipped_fib['N'] * np.log2(clipped_fib['N'])), 'b-', ls='dashed', label="Fibonacci L3 misses")
#ax2.plot(clipped_fib['N'], clipped_fib['instructions'] / (clipped_fib['N'] * np.log2(clipped_fib['N'])), 'b-', ls='dashed', label="Fib instructions")

plt.legend()
plt.savefig('runningtime.pdf')


# Instruction / cache miss plots
plt.figure()
ax1 = plt.subplot(111)
plt.title("Instructions and cache misses")
plt.xscale('log', basex=2)
ax1.set_ylabel('Normalized instruction count')
plt.xlabel('Length')

ax2 = ax1.twinx()
ax2.set_ylabel('Normalized L2 cache misses')

ax1.plot(clipped_fib['N'], clipped_fib['instructions'] / (clipped_fib['N'] * np.log2(clipped_fib['N'])), label="Fibonacci")
ax1.plot(clipped_hg_repetitive['N'], clipped_hg_repetitive['instructions'] / (clipped_hg_repetitive['N'] * np.log2(clipped_hg_repetitive['N'])), label="HG Repetitive")
ax1.plot(clipped_random['N'], clipped_random['instructions'] / (clipped_random['N'] * np.log2(clipped_random['N'])), label="Random")

ax2.plot(clipped_fib['N'], clipped_fib['L2_miss'] / (clipped_fib['N'] * np.log2(clipped_fib['N'])), ls='dashed', label="Fibonacci")
ax2.plot(clipped_hg_repetitive['N'], clipped_hg_repetitive['L2_miss'] / (clipped_hg_repetitive['N'] * np.log2(clipped_hg_repetitive['N'])), ls='dashed', label="HG Repetitive")
ax2.plot(clipped_random['N'], clipped_random['L2_miss'] / (clipped_random['N'] * np.log2(clipped_random['N'])), ls='dashed', label="HG Repetitive")

ax1.legend()
plt.savefig('instructions.pdf')
