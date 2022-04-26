import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import os
import json

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("file_path", type=Path)

p = parser.parse_args()

assert p.file_path.exists() , "The given directory doesn't exist"

output_dir = str(p.file_path)

if   output_dir.find('Example1a') != -1:
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example1a'
   ex_type = 1
elif output_dir.find('Example1b') != -1:
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example1b'
   ex_type = 2
elif output_dir.find('Example1c') != -1:
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example1c'
   ex_type = 3
elif output_dir.find('Example2c') != -1:
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example2c'
   ex_type = 4
elif output_dir.find('Example3a') != -1:
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example3a'
   ex_type = 5
else:
   print("Provided output directory: ", output_dir)
   print("*** Error: output directory path not a recognized format, quitting.")
   exit()

with open(output_dir+'/prob_config.json') as f:
   prob_config = json.load(f)


opt_type = prob_config['optimization_type'] 

if opt_type == 'dof_threshold':
   minimum_budget_problem = False
else:
   minimum_budget_problem = True


sns.set()
sns.set_context("talk", font_scale=3)
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')
palette_list = sns.color_palette(palette="tab10", n_colors=10)

def letterbox_entry(legend):
    from matplotlib.patches import Patch
    ax = legend.axes

    handles, labels = ax.get_legend_handles_labels()
    handles.insert(0, Patch(facecolor=palette_list[3], edgecolor='w'))
    labels.insert(0, "AMR policy costs")

    legend._legend_box = None
    legend._init_legend_box(handles, labels)
    legend._set_loc(legend._loc)
    legend.set_title(legend.get_title().get_text())


save_figs = True
# minimum_budget_problem = False

print("*** Check that correct data was loaded here in plots.py ***")
train_data_file = output_dir+'/training_data.csv'
rldata_file = output_dir+'/rl_data.csv'
deterministic_amr_data_file = output_dir+'/deterministic_amr_data.csv'
# rldata_file = output_dir+'/rl_data_ood.csv'
# deterministic_amr_data_file = output_dir+'/deterministic_amr_data_ood.csv'


df = pd.read_csv(train_data_file, index_col=0)
cost = df.index.to_numpy()

if ex_type == 4: # example 2c
   df = pd.read_csv(rldata_file)
   rlepisode_cost = df['rlepisode_cost'][0]
   rlactions = df[{'theta', 'rho'}].to_numpy()
   rlactions = rlactions[:-1,:] # drop last row, which was padded
   rldofs =  df['rldofs'].to_numpy()
   rlerrors =  df['rlerrors'].to_numpy()
else:
   df = pd.read_csv(rldata_file, index_col=0)
   rlactions = df.index.to_numpy()
   rlepisode_cost = rlactions[-1]
   rlactions = rlactions[:-1]
   rldofs = df.iloc[:, 0].to_numpy()
   rlerrors = df.iloc[:, 1].to_numpy()

df = pd.read_csv(deterministic_amr_data_file, index_col=0)

actions = df.index.to_numpy()
costs = df.iloc[:, 0].to_numpy()
errors = df.iloc[:, 1].to_numpy()
dofs = df.iloc[:, 2].to_numpy()

for i in range(len(dofs)):
   dofs[i] = np.array(eval(dofs[i]))
   errors[i] = np.array(eval(errors[i]))

plt.figure(figsize=(6,6))
plt.plot(cost[:-1], color=palette_list[4], lw=2, label=r'Training curve')
ax1 = plt.gca()
ax1.set_xlabel("Epoch")
if minimum_budget_problem:
   # ax1.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
   ax1.set_ylabel(r'mean episode cost')
else:
   # ax1.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
   ax1.set_ylabel(r'mean episode cost')
# ax1.legend()
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_training-curve.pdf',format='pdf', bbox_inches='tight')

# ## Make letter-box plot
plt.figure(figsize=(6,6))
ax2 = sns.boxenplot(y=costs, width=.6, color=palette_list[3], label='_nolegend_')
x2 = ax2.get_xlim()
plt.hlines(rlepisode_cost, x2[0], x2[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy cost')
y2 = ax2.get_ylim()
# plt.fill_between(x2, np.floor(y2[0]), rlepisode_cost, color=palette_list[9], label=r'Apparent performance barrier')
plt.fill_between(x2, np.floor(y2[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')

ax2.set_xlabel(r'')
if minimum_budget_problem:
   # ax2.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
   ax2.set_ylabel(r'cost at final step')
else:
   # ax2.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
   ax2.set_ylabel(r'cost at final step')
lgd = ax2.legend()
letterbox_entry(lgd)
sns.despine(ax=ax2, bottom=True)
ax2.tick_params(bottom=False)
plt.tight_layout()
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_rlepisode_cost.pdf',format='pdf', bbox_inches='tight')

## Plot theta vs. cost
plt.figure(figsize=(6,6))
ax3 = plt.gca()
plt.plot(actions[9::10], costs[9::10], 'o', color=palette_list[3], label=r'AMR policies')
x3 = ax3.get_xlim()
plt.hlines(rlepisode_cost, x3[0], x3[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy')
y3 = ax3.get_ylim()
# plt.fill_between(x3, np.floor(y3[0]), rlepisode_cost, color=palette_list[9], label=r'Apparent performance barrier')
plt.fill_between(x3, np.floor(y3[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')


ax3.set_xlabel(r'$\theta$ (constant) in AMR policy')
if minimum_budget_problem:
   # ax3.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
   ax3.set_ylabel(r'cost at final step')
else:
   # ax3.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
   ax3.set_ylabel(r'cost at final step')
ax3.legend(loc='upper center')
plt.tight_layout()

if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_fig3.pdf',format='pdf', bbox_inches='tight')

# ## Make convergence plots (1/2)
plt.figure(figsize=(6,6))
ax4 = plt.gca()
alpha = 1.0
plt.loglog(dofs[9],errors[9],'-o',lw=1.3, color=palette_list[3], alpha=alpha, label=r'AMR policies')
# plt.loglog(dofs[9][-1],errors[9][-1], marker="o", markersize=10, color=palette_list[3], label='_nolegend_')
for k in range(19,len(errors),10):
   plt.loglog(dofs[k],errors[k],'-o',lw=1.3, color=palette_list[3], label='_nolegend_')
   plt.loglog(dofs[k][-1],errors[k][-1], marker="o", markersize=10, color=palette_list[3], alpha=alpha, label='_nolegend_')
plt.loglog(rldofs,rlerrors,'-o',lw=1.3, color=palette_list[0], label=r'(AM)$^2$R policy')
plt.loglog(rldofs[-1],rlerrors[-1], marker="o", markersize=10, color=palette_list[0], label='_nolegend_')
ax4.set_xlabel(r'Degrees of freedom')
ax4.set_ylabel(r'Relative error')
ax4.legend()
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_fig4.pdf',format='pdf', bbox_inches='tight')

## Make convergence plots (2/2)
cumdofs = []
cumrldofs = np.cumsum(rldofs)
for k in range(len(dofs)):
   cumdofs.append(np.cumsum(dofs[k]))
plt.figure(figsize=(6,6))
ax5 = plt.gca()
plt.loglog(cumdofs[9],errors[9],'-o',lw=1.3, color=palette_list[3], alpha=alpha, label=r'AMR policies')
# plt.loglog(cumdofs[9][-1],errors[9][-1], marker="o", markersize=10, color=palette_list[3], label='_nolegend_')
for k in range(19,len(errors),10):
   plt.loglog(cumdofs[k],errors[k],'-o',lw=1.3, color=palette_list[3], label='_nolegend_')
   # plt.loglog(cumdofs[k][-1],errors[k][-1], marker="o", markersize=10, color=palette_list[3], alpha=alpha, label='_nolegend_')
plt.loglog(cumrldofs,rlerrors,'-o',lw=1.3, color=palette_list[0], label=r'(AM)$^2$R policy')
# plt.loglog(cumrldofs[-1],rlerrors[-1], marker="o", markersize=10, color=palette_list[0], label='_nolegend_')
ax5.set_xlabel(r'Cumulative degrees of freedom')
ax5.set_ylabel(r'Relative error')
ax5.legend()
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_fig5.pdf',format='pdf', bbox_inches='tight')

## Plot action vs. refinement step
plt.figure(figsize=(6,6))
ax6 = plt.gca()
plt.plot(rlactions,'-o',lw=1.3, label=r'(AM)$^2$R policy')
ax6.set_xlabel(r'Refinement step')
ax6.set_ylabel(r'$\theta$ in trained (AM)$^2$R policy')
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_fig6.pdf',format='pdf', bbox_inches='tight')


if ex_type == 4: # example 2c

   tp_data_file = output_dir+'/two_param_amr_data.csv'
   df = pd.read_csv(tp_data_file)

   plt.figure(figsize=(6,6))
   ax7 = sns.boxenplot(y=df['costs'], width=.6, color=palette_list[3], label='_nolegend_')
   x7 = ax7.get_xlim()
   plt.hlines(rlepisode_cost, x2[0], x2[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy cost')
   y2 = ax7.get_ylim()
   # plt.fill_between(x2, np.floor(y2[0]), rlepisode_cost, color=palette_list[9], label=r'Apparent performance barrier')
   plt.fill_between(x2, np.floor(y2[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')

   ax7.set_xlabel(r'')
   if minimum_budget_problem:
      # ax7.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
      ax7.set_ylabel(r'cost at final step')
   else:
      # ax7.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
      ax7.set_ylabel(r'cost at final step')
   lgd = ax7.legend()
   letterbox_entry(lgd)
   sns.despine(ax=ax7, bottom=True)
   ax7.tick_params(bottom=False)
   plt.tight_layout()
   if save_figs:
      plt.savefig(output_dir+'/'+fig_name_prefix+'_rl_vs_two_param.pdf',format='pdf', bbox_inches='tight')


plt.show()