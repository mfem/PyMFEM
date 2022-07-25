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
parser.add_argument("--second_file", type=Path)
parser.add_argument("--save_figs", default=True, action = 'store_true')
parser.add_argument('--no_save_figs', dest='savefigs', action='store_false')
parser.add_argument('--colorbar', default= True, action = 'store_true')
parser.add_argument('--no_colorbar', dest='colorbar', action='store_false') # label each alpha instead of using colorbar

args = parser.parse_args()
colorbar = args.colorbar
assert args.file_path.exists() , "The given file doesn't exist"
input_file = str(args.file_path)

# Check for second file to be plotted
try:
    assert args.second_file.exists(), "The second file doesn't exist"
    second_file = True
    input_file2 = str(args.second_file)
    
    # only allow second file with colorbar option
    assert colorbar == True, "Two files may only be plotted with the colorbar option"
except:
    second_file = False


# check that const_alpha_data.txt is the input file
assert input_file.find("const_alpha_data") != -1, "Error, no file const_alpha_data"
if second_file:
    assert input_file2.find("const_alpha_data") != -1, "Error, no second file const_alpha_data"


# determine experiment type, if any
if input_file.find("exp2.") != -1:
    exp_flag = 2
    label1 = 'no budget'
elif input_file.find("exp4.") != -1:
    exp_flag = 4
    label1 = 'with budget'
else:
    exp_flag = 0
    label1 = str(exp_flag)

file = open(input_file, "r")
lines = file.readlines()
n = len(lines) - 1

if second_file:
    file2 = open(input_file2, "r")
    lines2 = file2.readlines()
    n2 = len(lines2) - 1

    if input_file2.find("exp2.") != -1:
        exp_flag2 = 2
        label2 = 'no budget'
    elif input_file2.find("exp4.") != -1:
        exp_flag2 = 4
        label2 = 'with budget'
    else:
        exp_flag2 = 0
        label2 = str(exp_flag2)

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

if colorbar == False:
	# determine spacing for labels
	y_spacing     = 1.08*np.ones(n)
	x_spacing     = np.ones(n)

	if exp_flag == 2:
	    y_spacing[-2] = 1
	    y_spacing[-1] = 0.92
	    y_spacing[1]  = 1.12
	    y_spacing[3]  = 1.12

	    x_spacing[n-3:n] = 1.1

	elif exp_flag == 4:
	    x_spacing[1] = 1.1
	    y_spacing[1] = 1

	    y_spacing[4] = 1.12

	    x_spacing[10] = 1.1
	    y_spacing[10] = 1.1

# plot data for each alpha
alpha = np.zeros(n); cum_dofs = np.zeros(n); error = np.zeros(n);
for i in range(1, n+1):
    data     = lines[i].split(', ')

    alpha   [i-1] =       data[0]
    cum_dofs[i-1] = float(data[1])
    error   [i-1] = float(data[2])
    

    if colorbar == False:
    	plt.loglog(cum_dofs, error, '.k')
    	plt.annotate(r"$\alpha$ = " + alpha[i-1], (cum_dofs[i-1]*x_spacing[i-1], error[i-1]*y_spacing[i-1]))

if second_file == True:
    alpha2 = np.zeros(n2); cum_dofs2 = np.zeros(n2); error2 = np.zeros(n2);

    for i in range(1, n2+1):
        data2 = lines2[i].split(', ')

        alpha2   [i-1] =       data2[0]
        cum_dofs2[i-1] = float(data2[1])
        error2   [i-1] = float(data2[2])

    plt.scatter(cum_dofs2, error2, c = alpha2, marker = 's', label = label2)

if colorbar == True:
	plt.scatter(cum_dofs, error, c = alpha, label = label1)
	plt.colorbar().set_label(label = r'$\alpha$',size=20,weight='bold')
    
plt.xscale('log')
plt.yscale('log')

# plot title
plt.title(r'Example 1a with Multi-Objective Cost and Fixed $\alpha$', fontdict = {'fontsize':24})
if exp_flag == 4:
    plt.title(r'Example 1a with MO Cost, Fixed $\alpha$, and Observed Budget', fontdict = {'fontsize':24})

ax5.set_xlabel(r'Cumulative degrees of freedom $(J_k)$', fontsize=22)
ax5.set_ylabel(r'Global error estimate $(\eta_k)$', fontsize=22)
ax5.tick_params(axis='x', which = 'major', labelsize=22)
ax5.tick_params(axis='x', which = 'minor', labelsize=12)
ax5.tick_params(axis='y', labelsize=22)
# ax5.legend(loc='upper right', prop={'size': 14})
ax5.set_yticks([10**-2, 10**-3])

if second_file == True:
    plt.legend()

if args.save_figs:
    name = input_file.split('.')[0]
    if second_file:
        name2 = input_file2.split('.')[0].split('const_alpha_data')[1]
        plot_name = 'plot_' + name + name2 + '.png'  
    else:
        plot_name = 'plot_' + name + '.png'
    
    plt.savefig(plot_name, format='png', bbox_inches='tight')
