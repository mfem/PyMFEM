'''
Plot data corresponding to MO cost function for each fixed value of alpha.
'''

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("file_path", type=Path)
parser.add_argument("--second_file", type=Path)
parser.add_argument("--fixed_theta_data", type = Path)
parser.add_argument("--save_figs", default=True, action = 'store_true')
parser.add_argument('--no_save_figs', dest='savefigs', action='store_false')
parser.add_argument('--parameter_sweep', default = False, action = 'store_true')


args = parser.parse_args()
assert args.file_path.exists() , "The given file doesn't exist"
input_file = str(args.file_path)

# Check for second file to be plotted
try:
    assert args.second_file.exists(), "The second file doesn't exist"
    second_file = True
    input_file2 = str(args.second_file)

except:
    second_file = False

# Check for fixed_theta_data file
try:
    assert args.fixed_theta_data.exists(), "The fixed-theta data file doesn't exist."
    theta_file = True
    fixed_theta_file = str(args.fixed_theta_data)
except:
    theta_file = False

file = open(input_file, "r")
lines = file.readlines()
n = len(lines)

if second_file:
    file2 = open(input_file2, "r")
    lines2 = file2.readlines()
    n2 = len(lines2) - 1


# Process data in fixed-theta file, if it is given
if theta_file:
    df = pd.read_csv(args.fixed_theta_data)
    targets    = df['target'].unique()
    min_dofs   = [];
    min_errors = [];
    min_thetas = [];


    for target in targets:
        fdf = df[df['target'] == target]
        index = fdf['dofs'].argmin()
        min_dofs.append(fdf.iloc[index]['dofs'])
        #min_errors.append(fdf.iloc[index]['error']) # would like to add this if we recorded final error data
        min_thetas.append(fdf.iloc[index]['theta'])

# create figure
sns.set()
sns.set_context("talk", font_scale=3)
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')
# palette_list = sns.color_palette(palette="tab10", n_colors=10)

plt.figure(figsize=(6,6))
ax5 = plt.gca()

# plot data for each tau
taus = np.zeros(n); cum_dofs = np.zeros(n); error = np.zeros(n);
for i in range(0, n):
    data     = lines[i].split(', ')

    taus    [i] = 2**float(data[0])
    cum_dofs[i] =    float(data[1])
    error   [i] =    float(data[2])

plt.scatter(cum_dofs, error, label = "RL policy")

if second_file == True:
    taus2 = np.zeros(n2); cum_dofs2 = np.zeros(n2); error2 = np.zeros(n2);

    for i in range(1, n2+1):
        data2 = lines2[i].split(', ')

        taus2    [i-1] =    2**data2[0]
        cum_dofs2[i-1] = float(data2[1])
        error2   [i-1] = float(data2[2])

    plt.scatter(cum_dofs2, error2, c = taus2, marker = 's', label = label2)
    plt.legend()

# plot results from best fixed-theta policies
if theta_file:
    plt.scatter(min_dofs, targets, label = "Best fixed-theta policies")
    plt.legend()

'''
# plot results for fixed theta cases for a range of thetas
if args.parameter_sweep:
    ps_file_name = 'fixed_theta_1a_AMR_data.txt'
    #assert ps_file_name.exists() , "The parameter sweep file does not exist."

    # read in data
    file = open(ps_file_name, "r")
    lines = file.readlines()
    theta    = np.zeros(9)
    cum_dofs = np.zeros(9)
    error    = np.zeros(9)
    for i in range(1,10):
        data = lines[i].split(', ')

        theta   [i-1] = data[1]
        cum_dofs[i-1] = data[2]
        error   [i-1] = data[3]
    plt.scatter(cum_dofs, error, c=theta, marker = '*', label='Fixed-theta results')
    file.close()


if args.parameter_sweep == False:
    plt.colorbar().set_label(label = r'$\tau$',size=20,weight='bold')
else:
    plt.colorbar().set_label(label = r'$\tau$ or $\theta$',size=20,weight='bold')

if second_file or args.parameter_sweep:
    plt.legend()
'''
plt.xscale('log')
plt.yscale('log')

# plot title
plt.title(r'RL Policy vs Best Fixed-Theta Policies', fontdict = {'fontsize':24})

ax5.set_xlabel(r'Cumulative degrees of freedom $(J_k)$', fontsize=22)
ax5.set_ylabel(r'Global error estimate $(\eta_k)$', fontsize=22)
ax5.tick_params(axis='x', which = 'major', labelsize=22)
ax5.tick_params(axis='x', which = 'minor', labelsize=12)
ax5.tick_params(axis='y', labelsize=22)
# ax5.legend(loc='upper right', prop={'size': 14})
#ax5.set_yticks([10**-2, 10**-3])


if args.save_figs:
    directory = '/'.join(input_file.split('/')[:-1]) + '/'
    name = input_file.split('/')[-1].split('.')[0]
    
    plot_name = 'plot_' + name + '.png'
    
    plt.savefig(directory + plot_name, format='png', bbox_inches='tight')
    print("Saved pareto front figure to {}".format(directory + plot_name))
