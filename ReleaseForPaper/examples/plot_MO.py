'''
Plot data corresponding to MO cost function for each fixed value of alpha.
'''

import numpy as np
from matplotlib import pyplot as plt

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("file_path", type=Path)
parser.add_argument("--save_figs", defaul=True, action = 'store_true')
parser.add_argument('--no-savefigs', dest='savefigs', action='store_false')

args = parser.parse_args()
assert args.file_path.exists() , "The given file doesn't exist"
input_file = str(args.file_path)


# check that const_alpha_data.txt is the input file
assert input_file.find("const_alpha_data.txt") != -1, "Error, no file const_alpha_data.txt"

file = open(input_file, "r")
lines = file.readlines()
n = len(lines) - 1

# create figure
plt.figure(figsize=(6,6))
ax5 = plt.gca()

# plot data for each alpha
for i in range(1, n+1):
	data     = lines[i].split(', ')

	alpha    =       data[0]
	cum_dofs = float(data[1])
	error    = float(data[2])

	plt.loglog(cum_dofs, error, marker="^", lw=1.3, color=palette_list[0], label = "alpha = " + alpha)

ax5.set_xlabel(r'Cumulative degrees of freedom $(J_k)$', fontsize=22)
ax5.set_ylabel(r'Global error estimate $(\eta_k)$', fontsize=22)
ax5.tick_params(axis='x', labelsize=22)
ax5.tick_params(axis='y', labelsize=22)
ax5.legend(loc='upper right', prop={'size': 14})

if save_figs:
  plt.savefig('plot_const_alpha.pdf', format='pdf', bbox_inches='tight')
