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
hg_repetitive_dist = pd.read_csv("data/lcs_blowup_fasta_hg_repetitive_fasta_hg_repetitive_DIST.dat", delim_whitespace=True)
hg_nonrepetitive_dist = pd.read_csv("data/lcs_blowup_fasta_hg_nonrepetitive_fasta_hg_nonrepetitive_DIST.dat", delim_whitespace=True)
hg_combined_dist = pd.read_csv("data/lcs_blowup_fasta_hg_combined_fasta_hg_combined_DIST.dat", delim_whitespace=True)
random_dist = pd.read_csv("data/lcs_blowup_random_random_DIST.dat", delim_whitespace=True)

plt.plot(fibonacci_dist['A_len'] + fibonacci_dist['B_len'], normalize_dist(fibonacci_dist), label="Fibonacci")
plt.plot(hg_repetitive_dist['A_len'] + hg_repetitive_dist['B_len'], normalize_dist(hg_repetitive_dist), label="HG Repetitive")
plt.plot(hg_nonrepetitive_dist['A_len'] + hg_nonrepetitive_dist['B_len'], normalize_dist(hg_nonrepetitive_dist), label="HG Non-repetitive")
plt.plot(hg_combined_dist['A_len'] + hg_combined_dist['B_len'], normalize_dist(hg_combined_dist), label="HG Combined")
plt.plot(random_dist['A_len'] + random_dist['B_len'], normalize_dist(random_dist), color='k', label="Random")

plt.legend()
plt.savefig('dist_runningtime.pdf')

# Instruction count and cache misses for DIST repository computation
plt.figure()
plt.title("Running time of DIST Repository computation")
ax1 = plt.subplot(111)
plt.xscale('log', basex=2)
plt.ylabel('Normalized instruction count')
plt.xlabel('N')

ax2 = ax1.twinx()
ax2.set_ylabel('L2 cache miss percentage')

# ax1.plot(fibonacci_dist['A_len'] + fibonacci_dist['B_len'], normalize_dist(fibonacci_dist, 'instructions'), label="Fibonacci")
ax1.plot(random_dist['A_len'] + random_dist['B_len'], normalize_dist(random_dist, 'instructions'), color='k', label="Random")
ax1.plot(hg_repetitive_dist['A_len'] + hg_repetitive_dist['B_len'], normalize_dist(hg_repetitive_dist, 'instructions'), color='r', label="HG Repetitive")
ax1.plot(hg_nonrepetitive_dist['A_len'] + hg_nonrepetitive_dist['B_len'], normalize_dist(hg_nonrepetitive_dist, 'instructions'), color='c', label="HG Non-repetitive")
ax1.plot(hg_combined_dist['A_len'] + hg_combined_dist['B_len'], normalize_dist(hg_combined_dist, 'instructions'), color='m', label="HG Non-repetitive")

ax2.plot(random_dist['A_len'] + random_dist['B_len'], (random_dist['L2_miss'] * 100) / (random_dist['L2_hits'] + random_dist['L2_miss']), label="Fibonacci", ls='dashed', color='k')
ax2.plot(hg_repetitive_dist['A_len'] + hg_repetitive_dist['B_len'], (hg_repetitive_dist['L2_miss'] * 100) / (hg_repetitive_dist['L2_hits'] + hg_repetitive_dist['L2_miss']), label="HG Repetitive", color='r', ls='dashed')
ax2.plot(hg_nonrepetitive_dist['A_len'] + hg_nonrepetitive_dist['B_len'], (hg_nonrepetitive_dist['L2_miss'] * 100) / (hg_nonrepetitive_dist['L2_hits'] + hg_nonrepetitive_dist['L2_miss']), label="HG Non-repetitive", color='c', ls='dashed')
ax2.plot(hg_combined_dist['A_len'] + hg_combined_dist['B_len'], (hg_combined_dist['L2_miss'] * 100) / (hg_combined_dist['L2_hits'] + hg_combined_dist['L2_miss']), label="HG Combined", color='m', ls='dashed')

ax1.legend()
plt.savefig('dist_runningtime_cpu.pdf')

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
plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,2))

fibonacci_grid = pd.read_csv("data/lcs_blowup_fib_fib_grid.dat", delim_whitespace=True)
# genome_grid = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_grid.dat", delim_whitespace=True)
hg_repetitive_grid = pd.read_csv("data/lcs_blowup_fasta_hg_repetitive_fasta_hg_repetitive_grid.dat", delim_whitespace=True)
hg_nonrepetitive_grid = pd.read_csv("data/lcs_blowup_fasta_hg_nonrepetitive_fasta_hg_nonrepetitive_grid.dat", delim_whitespace=True)
hg_combined_grid = pd.read_csv("data/lcs_blowup_fasta_hg_combined_fasta_hg_combined_grid.dat", delim_whitespace=True)
genome_grid = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_grid.dat", delim_whitespace=True)
random_grid = pd.read_csv("data/lcs_blowup_random_random_grid.dat", delim_whitespace=True)

plt.plot(fibonacci_grid['A_len'] + fibonacci_grid['B_len'], normalize_grid(fibonacci_grid), label="Fibonacci")
# plt.plot(genome_grid['A_len'] + genome_grid['B_len'], normalize_grid(genome_grid), label="Genome 1 - 2")
plt.plot(hg_repetitive_grid['A_len'] + hg_repetitive_grid['B_len'], normalize_grid(hg_repetitive_grid), label="HG Repetitive")
plt.plot(hg_nonrepetitive_grid['A_len'] + hg_nonrepetitive_grid['B_len'], normalize_grid(hg_nonrepetitive_grid), label="HG Nonrepetitive")
plt.plot(hg_combined_grid['A_len'] + hg_combined_grid['B_len'], normalize_grid(hg_combined_grid), label="HG Combined")
plt.plot(random_grid['A_len'] + random_grid['B_len'], normalize_grid(random_grid), color='k', label="Random")

plt.legend()
plt.savefig('grid_runningtime.pdf')


# CPU stats for grid computation
plt.figure()
plt.title("CPU statistics for the grid computation")
plt.xscale('log', basex=2)
plt.xlabel('N')
plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,2))

ax1 = plt.subplot(111)
ax1.set_ylabel('Normalized instruction count')
ax2 = ax1.twinx()
ax2.set_ylabel('Normalized L3 cache misses')

ax1.plot((fibonacci_grid['A_len'] + fibonacci_grid['B_len']).tail(23), normalize_grid(fibonacci_grid, 'instructions').tail(23), color='b', label="Fibonacci")
ax2.plot((fibonacci_grid['A_len'] + fibonacci_grid['B_len']).tail(23), normalize_grid(fibonacci_grid, 'L3_miss').tail(23), ls='dashed', color='b', label="Fibonacci")

ax1.legend(loc='upper center')
plt.savefig('grid_runningtime_cpu.pdf')


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
hg_repetitive_total = pd.read_csv("data/lcs_blowup_fasta_hg_repetitive_fasta_hg_repetitive_total.dat", delim_whitespace=True)
hg_nonrepetitive_total = pd.read_csv("data/lcs_blowup_fasta_hg_nonrepetitive_fasta_hg_nonrepetitive_total.dat", delim_whitespace=True)
hg_combined_total = pd.read_csv("data/lcs_blowup_fasta_hg_combined_fasta_hg_combined_total.dat", delim_whitespace=True)
genome_total = pd.read_csv("data/lcs_blowup_fasta_genome1_fasta_genome2_total.dat", delim_whitespace=True)

plt.plot(fibonacci_total['A_len'] + fibonacci_total['B_len'], normalize_total(fibonacci_total), label="Fibonacci")
plt.plot(hg_repetitive_total['A_len'] + hg_repetitive_total['B_len'], normalize_total(hg_repetitive_total), label="HG Repetitive")
plt.plot(hg_nonrepetitive_total['A_len'] + hg_nonrepetitive_total['B_len'], normalize_total(hg_nonrepetitive_total), label="HG Non-repetitive")
plt.plot(hg_combined_total['A_len'] + hg_combined_total['B_len'], normalize_total(hg_combined_total), ls='dashed', label="HG Combined")
plt.plot(random_total['A_len'] + random_total['B_len'], normalize_total(random_total), color='k', label="Random")

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


# Plot of the relative time consumption between the different parts of the algorithm for hg_repetitive input
fig = plt.figure()
plt.title("Time consumed by different steps of the algorithm")
plt.xscale('log', basex=2)
plt.ylabel('Relative time used')
plt.xlabel('N')

hg_repetitive_slp = pd.read_csv("data/lcs_blowup_fasta_hg_repetitive_fasta_hg_repetitive_SLP.dat", delim_whitespace=True)
N = (hg_repetitive_total['A_len'] + hg_repetitive_total['B_len']).values
slp_time = ((hg_repetitive_slp['median'] / hg_repetitive_total['median']) * 100).values
dist_time = ((hg_repetitive_dist['median'] / hg_repetitive_total['median']) * 100).values
grid_time = ((hg_repetitive_grid['median'] / hg_repetitive_total['median']) * 100).values

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

plt.savefig('hg_repetitive_area_plot.pdf')


# Plot of the relative time consumption between the different parts of the algorithm for hg_nonrepetitive input
fig = plt.figure()
plt.title("Time consumed by different steps of the algorithm")
plt.xscale('log', basex=2)
plt.ylabel('Relative time used')
plt.xlabel('N')

hg_nonrepetitive_slp = pd.read_csv("data/lcs_blowup_fasta_hg_nonrepetitive_fasta_hg_nonrepetitive_SLP.dat", delim_whitespace=True)
N = (hg_nonrepetitive_total['A_len'] + hg_nonrepetitive_total['B_len']).values
slp_time = ((hg_nonrepetitive_slp['median'] / hg_nonrepetitive_total['median']) * 100).values
dist_time = ((hg_nonrepetitive_dist['median'] / hg_nonrepetitive_total['median']) * 100).values
grid_time = ((hg_nonrepetitive_grid['median'] / hg_nonrepetitive_total['median']) * 100).values

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

plt.savefig('hg_nonrepetitive_area_plot.pdf')


# Plot of the relative time consumption between the different parts of the algorithm for hg_combined input
fig = plt.figure()
plt.title("Time consumed by different steps of the algorithm")
plt.xscale('log', basex=2)
plt.ylabel('Relative time used')
plt.xlabel('N')

hg_combined_slp = pd.read_csv("data/lcs_blowup_fasta_hg_combined_fasta_hg_combined_SLP.dat", delim_whitespace=True)
N = (hg_combined_total['A_len'] + hg_combined_total['B_len']).values
slp_time = ((hg_combined_slp['median'] / hg_combined_total['median']) * 100).values
dist_time = ((hg_combined_dist['median'] / hg_combined_total['median']) * 100).values
grid_time = ((hg_combined_grid['median'] / hg_combined_total['median']) * 100).values

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

plt.savefig('hg_combined_area_plot.pdf')


# Running time of Simple vs LCSBlowup
simple_fibonacci = pd.read_csv("data/simple_fib_fib_total.dat", delim_whitespace=True)
simple_genome = pd.read_csv("data/simple_fasta_genome1_fasta_genome2_total.dat", delim_whitespace=True)
simple_hg_repetitive = pd.read_csv("data/simple_fasta_hg_repetitive_fasta_hg_repetitive_total.dat", delim_whitespace=True)
simple_hg_nonrepetitive = pd.read_csv("data/simple_fasta_hg_nonrepetitive_fasta_hg_nonrepetitive_total.dat", delim_whitespace=True)
simple_hg_combined = pd.read_csv("data/simple_fasta_hg_combined_fasta_hg_combined_total.dat", delim_whitespace=True)
simple_random = pd.read_csv("data/simple_random_random_total.dat", delim_whitespace=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

plt.title('Running time of Simple and LCSBlowUp algorithms')
plt.xscale('log', basex=2)
plt.yscale('log', basey=2)
plt.xlabel('N')
plt.ylabel('Speedup factor')

plt.plot(fibonacci_total['A_len'] + fibonacci_total['B_len'], simple_fibonacci['median'] / fibonacci_total['median'], label="Fibonacci")
# plt.plot(genome_total['A_len'] + genome_total['B_len'], simple_genome['median'] / genome_total['median'], label="Genomre 1-2")
plt.plot(hg_repetitive_total['A_len'] + hg_repetitive_total['B_len'], simple_hg_repetitive['median'] / hg_repetitive_total['median'], label="HG Repetitive")
plt.plot(hg_nonrepetitive_total['A_len'] + hg_nonrepetitive_total['B_len'], simple_hg_nonrepetitive['median'] / hg_nonrepetitive_total['median'], label="HG Non-repetitive")
plt.plot(hg_combined_total['A_len'] + hg_combined_total['B_len'], simple_hg_combined['median'] / hg_combined_total['median'], label="HG Combined")
plt.plot(random_total['A_len'] + random_total['B_len'], simple_random['median'] / random_total['median'], color='k', label="Random")

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                  box.width, box.height * 0.85])

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 3)
plt.savefig('simple_vs_lcs.pdf')
