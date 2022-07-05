'''
Plot data corresponding to MO cost function for each fixed value of alpha.
'''

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("file_path", type=Path)
parser.add_argument("--save_figs", default=True, action = 'store_true')
parser.add_argument('--no_save_figs', dest='savefigs', action='store_false')

args = parser.parse_args()
assert args.file_path.exists() , "The given file doesn't exist"
input_file = str(args.file_path)


# check that const_alpha_data.txt is the input file
assert input_file.find("const_alpha_data.txt") != -1, "Error, no file const_alpha_data.txt"

file = open(input_file, "r")
lines = file.readlines()
n = len(lines) - 1

# create figure
sns.set()
sns.set_context("talk", font_scale=3)
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')
# palette_list = sns.color_palette(palette="tab10", n_colors=10)

plt.figure(figsize=(6,6))
ax5 = plt.gca()

# determine spacing for labels
y_spacing     = 1.08*np.ones(n)
y_spacing[-2] = 1
y_spacing[-1] = 0.92
y_spacing[1]  = 1.12
y_spacing[3]  = 1.12

x_spacing     = np.ones(n)
x_spacing[n-3:n] = 1.1

# plot data for each alpha
for i in range(1, n+1):
    data     = lines[i].split(', ')

    alpha    =       data[0]
    cum_dofs = float(data[1])
    error    = float(data[2])

    plt.loglog(cum_dofs, error, '.k')
    plt.annotate(r"$\alpha$ = " + alpha, (cum_dofs*x_spacing[i-1], error*y_spacing[i-1]))

plt.xscale('log')
plt.yscale('log')
plt.title(r'Example 1a with Multi-Objective Cost and Fixed $\alpha$', fontdict = {'fontsize':24})
ax5.set_xlabel(r'Cumulative degrees of freedom $(J_k)$', fontsize=22)
ax5.set_ylabel(r'Global error estimate $(\eta_k)$', fontsize=22)
ax5.tick_params(axis='x', which = 'major', labelsize=22)
ax5.tick_params(axis='x', which = 'minor', labelsize=15)
ax5.tick_params(axis='y', labelsize=22)
# ax5.legend(loc='upper right', prop={'size': 14})
ax5.set_yticks([10**-2, 10**-3])


if args.save_figs:
  plt.savefig('plot_const_alpha.pdf', format='pdf', bbox_inches='tight')
