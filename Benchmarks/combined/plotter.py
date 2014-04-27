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

# Plot to verify the theoretical running time of the DIST repository computation
def normalize_dist(df, to_normalize = 'median'):
  df_x = df['A_x'] + df['B_x']
  df_n = df['A_prod'] + df['B_prod']
  return df[to_normalize] / (np.power(df_n, 2) * (df_x * np.log2(df_x)))

plt.figure()
plt.title("Running time of DIST Repository computation")
plt.xscale('log', basex=2)
plt.ylabel('Normalized running time')
plt.xlabel('N')

fibonacci_dist = pd.read_csv("data/lcs_blowup_fib_fib_DIST.dat", delim_whitespace=True)
genome_dist = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_DIST.dat", delim_whitespace=True)
random_dist = pd.read_csv("data/lcs_blowup_random_random_DIST.dat", delim_whitespace=True)

plt.plot(fibonacci_dist['A_len'] + fibonacci_dist['B_len'], normalize_dist(fibonacci_dist), label="Fibonacci")
plt.plot(genome_dist['A_len'] + genome_dist['B_len'], normalize_dist(genome_dist), label="Genome 1 - 2")
plt.plot(random_dist['A_len'] + random_dist['B_len'], normalize_dist(random_dist), label="Random")

plt.legend()
plt.savefig('dist_runningtime.pdf')


# Plot to verify the theoretical running time for filling out the grid
def normalize_grid(df, to_normalize = 'median'):
  df_x = df['A_x'] + df['B_x']
  df_N = df['A_len'] + df['B_len']
  return df[to_normalize] / (np.power(df_N, 2) / df_x)

plt.figure()
plt.title("Running time for the grid computation")
plt.xscale('log', basex=2)
plt.ylabel('Normalized running time')
plt.xlabel('N')

fibonacci_grid = pd.read_csv("data/lcs_blowup_fib_fib_grid.dat", delim_whitespace=True)
genome_grid = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_grid.dat", delim_whitespace=True)
random_grid = pd.read_csv("data/lcs_blowup_random_random_grid.dat", delim_whitespace=True)

plt.plot(fibonacci_grid['A_len'] + fibonacci_grid['B_len'], normalize_grid(fibonacci_grid), label="Fibonacci")
plt.plot(genome_grid['A_len'] + genome_grid['B_len'], normalize_grid(genome_grid), label="Genome 1 - 2")
plt.plot(random_grid['A_len'] + random_grid['B_len'], normalize_grid(random_grid), label="Random")

plt.legend()
plt.savefig('grid_runningtime.pdf')


# Plot to verify theoretical running time of the combined algorithm
def normalize_total(df, to_normalize_total = 'median'):
  return df[to_normalize_total] / ((df['A_len'] + df['B_len']) * (df['A_prod'] + df['B_prod']) * np.log2(np.max((df['A_len'] + df['B_len']) / (df['A_prod'] + df['B_prod']), 2)))

plt.figure()
plt.title("Running time of combined algorithm")
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)
plt.ylabel('Normalized running time')
plt.xlabel('N')

fibonacci_total = pd.read_csv("data/lcs_blowup_fib_fib_total.dat", delim_whitespace=True)
random_total = pd.read_csv("data/lcs_blowup_random_random_total.dat", delim_whitespace=True)
genome_total = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_total.dat", delim_whitespace=True)

plt.plot(fibonacci_total['A_len'] + fibonacci_total['B_len'], normalize_total(fibonacci_total), label="Fibonacci")
plt.plot(genome_total['A_len'] + genome_total['B_len'], normalize_total(genome_total), label="Genome 1 - 2")
plt.plot(random_total['A_len'] + random_total['B_len'], normalize_total(random_total), label="Random")

plt.legend()
plt.savefig('total_runningtime.pdf')


# Plot of the relative time consumption between the different parts of the algorithm for Fibonacci input
fig = plt.figure()
plt.title("Time consumed by different steps of the algorithm")
plt.xscale('log', basex=2)
plt.ylabel('Relative time used')
plt.xlabel('N')

fibonacci_slp = pd.read_csv("data/lcs_blowup_fib_fib_SLP.dat", delim_whitespace=True)
N = (fibonacci_total['A_len'] + fibonacci_total['B_len']).values
slp_time = ((fibonacci_slp['median'] / fibonacci_total['median']) * 100).values
dist_time = ((fibonacci_dist['median'] / fibonacci_total['median']) * 100).values
grid_time = ((fibonacci_grid['median'] / fibonacci_total['median']) * 100).values

y = np.row_stack((slp_time, dist_time, grid_time))
y_stack = np.cumsum(y, axis=0)

ax1 = fig.add_subplot(111)
ax1.fill_between(N, 0, y_stack[0,:], facecolor="#AC4040")
ax1.fill_between(N, y_stack[0,:], y_stack[1,:], facecolor="#1DACD6")
ax1.fill_between(N, y_stack[1,:], y_stack[2,:], facecolor="#1DD6AC")

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])

proxy1 = plt.Rectangle((0, 0), 1, 1, facecolor="#AC4040")
proxy2 = plt.Rectangle((0, 0), 1, 1, facecolor="#1DACD6")
proxy3 = plt.Rectangle((0, 0), 1, 1, facecolor="#1DD6AC")
plt.legend([proxy1, proxy2, proxy3], ['SLP', 'DIST Repository', 'Fill out grid'], loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 3)

plt.savefig('fib_area_plot.pdf')


# Plot of the relative time consumption between the different parts of the algorithm for random input
fig = plt.figure()
plt.title("Time consumed by different steps of the algorithm")
plt.xscale('log', basex=2)
plt.ylabel('Relative time used')
plt.xlabel('N')

random_slp = pd.read_csv("data/lcs_blowup_random_random_SLP.dat", delim_whitespace=True)
N = (random_total['A_len'] + random_total['B_len']).values
slp_time = ((random_slp['median'] / random_total['median']) * 100).values
dist_time = ((random_dist['median'] / random_total['median']) * 100).values
grid_time = ((random_grid['median'] / random_total['median']) * 100).values

y = np.row_stack((slp_time, dist_time, grid_time))
y_stack = np.cumsum(y, axis=0)

ax1 = fig.add_subplot(111)
ax1.fill_between(N, 0, y_stack[0,:], facecolor="#AC4040")
ax1.fill_between(N, y_stack[0,:], y_stack[1,:], facecolor="#1DACD6")
ax1.fill_between(N, y_stack[1,:], y_stack[2,:], facecolor="#1DD6AC")

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])

proxy1 = plt.Rectangle((0, 0), 1, 1, facecolor="#AC4040")
proxy2 = plt.Rectangle((0, 0), 1, 1, facecolor="#1DACD6")
proxy3 = plt.Rectangle((0, 0), 1, 1, facecolor="#1DD6AC")
plt.legend([proxy1, proxy2, proxy3], ['SLP', 'DIST Repository', 'Fill out grid'], loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 3)

plt.savefig('random_area_plot.pdf')


# Plot of the relative time consumption between the different parts of the algorithm for genome input
fig = plt.figure()
plt.title("Time consumed by different steps of the algorithm")
plt.xscale('log', basex=2)
plt.ylabel('Relative time used')
plt.xlabel('N')

genome_slp = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_SLP.dat", delim_whitespace=True)
N = (genome_total['A_len'] + genome_total['B_len']).values
slp_time = ((genome_slp['median'] / genome_total['median']) * 100).values
dist_time = ((genome_dist['median'] / genome_total['median']) * 100).values
grid_time = ((genome_grid['median'] / genome_total['median']) * 100).values

y = np.row_stack((slp_time, dist_time, grid_time))
y_stack = np.cumsum(y, axis=0)

ax1 = fig.add_subplot(111)
ax1.fill_between(N, 0, y_stack[0,:], facecolor="#AC4040")
ax1.fill_between(N, y_stack[0,:], y_stack[1,:], facecolor="#1DACD6")
ax1.fill_between(N, y_stack[1,:], y_stack[2,:], facecolor="#1DD6AC")

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])

proxy1 = plt.Rectangle((0, 0), 1, 1, facecolor="#AC4040")
proxy2 = plt.Rectangle((0, 0), 1, 1, facecolor="#1DACD6")
proxy3 = plt.Rectangle((0, 0), 1, 1, facecolor="#1DD6AC")
plt.legend([proxy1, proxy2, proxy3], ['SLP', 'DIST Repository', 'Fill out grid'], loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 3)

plt.savefig('genome_area_plot.pdf')


# Running time of Simple vs LCSBlowup
simple_fibonacci = pd.read_csv("data/simple_fib_fib_total.dat", delim_whitespace=True)
simple_genome = pd.read_csv("data/simple_fasta_genome1_fasta_genome2_total.dat", delim_whitespace=True)
simple_random = pd.read_csv("data/simple_random_random_total.dat", delim_whitespace=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

plt.title('Running time of Simple and LCSBlowUp algorithms')
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)
plt.xlabel('N')
plt.ylabel('Speedup factor')

plt.plot(fibonacci_total['A_len'] + fibonacci_total['B_len'], simple_fibonacci['median'] / fibonacci_total['median'], label="Fibonacci")
plt.plot(genome_total['A_len'] + genome_total['B_len'], simple_genome['median'] / genome_total['median'], label="Genomre 1-2")
plt.plot(random_total['A_len'] + random_total['B_len'], simple_random['median'] / random_total['median'], label="Random")

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 3)
plt.savefig('simple_vs_lcs.pdf')
