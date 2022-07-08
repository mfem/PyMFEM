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
parser.add_argument('--mesh_abbrv', type=str, required=False)
parser.add_argument('--angle_abbrv', type=str, required=False)

args = parser.parse_args()

assert args.file_path.exists() , "The given directory doesn't exist"

output_dir = str(args.file_path)

if args.angle_abbrv is not None:
   print("Plotting for angle ", args.angle_abbrv)

if   output_dir.find('Example1a_MO') != -1:
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example1a_MO'
   ex_type = 6

   # determine if alpha was fixed 
   if output_dir.find('alpha') != -1:
       alpha_fixed = True
       
       # determine value of alpha from directory name 
       alpha_str = output_dir.split('alpha')[1].split('_')
       alpha_str = alpha_str[0] + '.' + alpha_str[1]

   else:
       alpha_fixed = False

elif output_dir.find('Example1a') != -1:
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
   assert args.angle_abbrv is not None , "Need to provide angle to plots.py for Example 2c"
   assert args.mesh_abbrv is not None, "Need to provide mesh name to plots.py for Example 2c: "
   print("Loading data from ", output_dir)
   fig_name_prefix = 'Example2c_' + args.mesh_abbrv + '_angle_' + args.angle_abbrv
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
have_expert_policy = False
if fig_name_prefix == 'Example1a' or fig_name_prefix == 'Example1a_MO':
   have_expert_policy = True

print("*** Check that correct data was loaded here in plots.py ***")
train_data_file = output_dir+'/training_data.csv'
rldata_file = output_dir+'/rl_data.csv'
expert_file = output_dir+'/expert_data.csv'
# expert_file = output_dir+'/deterministic_amr_data.csv'
if args.angle_abbrv is not None:
   rldata_file = rldata_file[:-4] + '_' + args.mesh_abbrv + '_angle_' + args.angle_abbrv + '.csv'
   expert_file = expert_file[:-4] + '_' + args.mesh_abbrv + '_angle_' + args.angle_abbrv + '.csv'
   twopar_file = output_dir + '/tpp_data_' + args.mesh_abbrv + '_angle_' + args.angle_abbrv + '.csv'
   print("RL data file = \n   ", rldata_file)
   print("two param data file = \n   ", twopar_file)
   if args.angle_abbrv == 'nan':
      print("expert data file = \n   ", expert_file)
      have_expert_policy = False



df = pd.read_csv(train_data_file, index_col=0)
cost = df.index.to_numpy()

if ex_type == 4: # example 2c
   df = pd.read_csv(rldata_file)
   rlepisode_cost = df['rlepisode_cost'][0]
   rlactions_theta = df['theta'].to_numpy()
   rlactions_rho   = df['rho'].to_numpy()
   rlactions_theta = rlactions_theta[:-1] # drop last entry, which was padded
   rlactions_rho = rlactions_rho[:-1] # drop last entry, which was padded
   rldofs =  df['rldofs'].to_numpy()
   rlerrors =  df['rlerrors'].to_numpy()
else:
   df = pd.read_csv(rldata_file, index_col=0)
   rlactions = df.index.to_numpy()
   rlepisode_cost = rlactions[-1]
   rlactions = rlactions[:-1] # drop last entry, which was padded
   rldofs = df.iloc[:, 0].to_numpy()
   rlerrors = df.iloc[:, 1].to_numpy()

if have_expert_policy:
   df = pd.read_csv(expert_file, index_col=0)
   actions = df.index.to_numpy()
   costs = df.iloc[:, 0].to_numpy()
   errors = df.iloc[:, 1].to_numpy()
   dofs = df.iloc[:, 2].to_numpy()
   for i in range(len(dofs)):
      dofs[i] = np.array(eval(dofs[i]))
      errors[i] = np.array(eval(errors[i]))


##########
# make fig training-curve
##########
plt.figure(figsize=(6,6))
plt.plot(cost[:-1], color=palette_list[4], lw=2, label=r'Training curve')
ax1 = plt.gca()
ax1.set_xlabel("Number of training batches", fontsize=22)
if minimum_budget_problem:
   # ax1.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
   ax1.set_ylabel(r'Mean episode cost $\left(\mathbb{E}_{\theta\sim\pi}\left[ \log_2(J_k)\right]\right)$', fontsize=22)
else:
   # ax1.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
   ax1.set_ylabel(r'Mean episode cost $\left(\mathbb{E}_{\theta\sim\pi}\left[ \log_2(\eta_k)\right]\right)$', fontsize=22)
# ax1.legend()
ax1.tick_params(axis='x', labelsize=22)
ax1.tick_params(axis='y', labelsize=22)
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_training-curve.pdf',format='pdf', bbox_inches='tight')


##########
# make fig6: Plot action vs. refinement step
##########

plt.figure(figsize=(6,6))
ax6 = plt.gca()
if ex_type == 4: # expand this "if" condition to accomdate any hp problem
   # print(rlactions[:,0])
   # plot(x1, y1, 'bo')
   plt.plot(rlactions_theta,'-o',ms=12.0, lw=4.0, ls='solid', color=palette_list[0], label=r'$\theta$ (h parameter)')
   plt.plot(rlactions_rho,  '-^',ms=12.0, lw=4.0, ls='solid', color=palette_list[1], label=r'$\rho$ (p parameter)')
   
   # # overlay best two parameter policy
   # df_tpp = pd.read_csv(twopar_file)
   # tpp_best_theta = df_tpp.iloc[df_tpp['costs'].argmin()]['theta'] * 1.002 * np.ones(len(rlactions_theta))
   # tpp_best_rho = df_tpp.iloc[df_tpp['costs'].argmin()]['rho'] * 0.998 * np.ones(len(rlactions_rho))
   # plt.plot(tpp_best_theta,'-x',lw=1.3, color=palette_list[0], label=r'optimal fixed $\theta$')
   # plt.plot(tpp_best_rho,'-x',lw=1.3, color=palette_list[1], label=r'optimal fixed $\rho$')

   # # overlay best two parameter policy according to meshes seen in training: theta = 0.6, rho = 0.4
   df_tpp = pd.read_csv(twopar_file)
   best_theta = df_tpp.iloc[df_tpp['costs'].argmin()]['theta']
   best_rho   = df_tpp.iloc[df_tpp['costs'].argmin()]['rho']
   tpp_ep_cost = df_tpp[df_tpp['theta'] == best_theta][df_tpp['rho'] == best_rho]['costs'].to_numpy().item()
   print("Best tpp ep cost = ", tpp_ep_cost, " occured at (theta, rho)=(",best_theta, ", ",best_rho,")")
   print("RL cost = ", rlepisode_cost)
   # plt.plot(0.6 * np.ones(len(rlactions_theta)),'-x',lw=1.3, color=palette_list[0], label=r'optimal fixed $\theta$')
   # plt.plot(0.5 * np.ones(len(rlactions_rho)),'-x',lw=1.3, color=palette_list[1], label=r'optimal fixed $\rho$')

   
   ax6.set_xlim(-0.5, 14.5) # forces x axis to have ticks from 0 to 14
   ax6.set_ylim(0.0, 0.8)
   plt.xticks(fontsize=30)
   plt.yticks(fontsize=30)
   policy_comp_string = 'mesh = ' + args.mesh_abbrv + ' angle = ' + args.angle_abbrv  + ' tpp cost = ' +  str(np.round(tpp_ep_cost,2)) + ' RL cost = ' + str(np.round(rlepisode_cost, 2))
   text_file = open("policy_comp.txt", "a")
   n = text_file.write(policy_comp_string + "\n")
   text_file.close()
   # ax6.set_xlabel('mesh = ' + args.mesh_abbrv + ' tpp cost = ' +  str(np.round(tpp_ep_cost,2)) + ' RL cost = ' + str(np.round(rlepisode_cost, 2)))
   # ax6.set_ylabel(r'parameter values in trained (AM)$^2$R policy')
   ax6.legend(loc='lower left', prop={'size': 26})

else:
   plt.plot(rlactions,'-o',lw=1.3, label=r'(AM)$^2$R policy')
   ax6.set_ylabel(r'$\theta$ selected by trained (AM)$^2$R policy', fontsize=22)
   ax6.set_xlabel(r'Refinement iteration $(k)$', fontsize=22)
   ax6.tick_params(axis='x', labelsize=22)
   ax6.tick_params(axis='y', labelsize=22)
if save_figs:
   plt.savefig(output_dir+'/'+fig_name_prefix+'_fig6.pdf',format='pdf', bbox_inches='tight')



if have_expert_policy:
   ##########
   # make fig3
   ##########
   ## Plot expert theta vs. cost
   plt.figure(figsize=(6,6))
   ax3 = plt.gca()
   plt.plot(actions[9::10], costs[9::10], 'o', color=palette_list[3], label=r'AMR policies')
   x3 = ax3.get_xlim()
   plt.hlines(rlepisode_cost, x3[0], x3[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy')
   y3 = ax3.get_ylim()
   # note: factor 0.1 below makes hatch area have height = 10% of the range of the y3 axis
   plt.fill_between(x3, y3[0] - 0.1 * (y3[1]-y3[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')
   ax3.set_xlabel(r'Fixed parameter in AMR policy ($\theta$)', fontsize=22)
   if minimum_budget_problem:
      ax3.set_ylabel(r'Cost at final step $(\log_2(J_k))$', fontsize=22)
   else:
      ax3.set_ylabel(r'Final error $\left(\log_2(\eta_k)\right)$', fontsize=22)
   ax3.tick_params(axis='x', labelsize=22)
   ax3.tick_params(axis='y', labelsize=22)
   ax3.legend(loc='upper left', prop={'size': 14})
   plt.tight_layout()

   if save_figs:
      plt.savefig(output_dir+'/'+fig_name_prefix+'_fig3.pdf',format='pdf', bbox_inches='tight')

   ##########
   # make fig4
   ##########
   # ## Make convergence plots (1/2)
   plt.figure(figsize=(6,6))
   ax4 = plt.gca()
   alpha = 1.0
   plt.loglog(dofs[9],errors[9],'-o',lw=1.3, color=palette_list[3], alpha=alpha, label=r'AMR policies')
   # plt.loglog(dofs[9][-1],errors[9][-1], marker="o", markersize=10, color=palette_list[3], label='_nolegend_')
   for k in range(19,len(errors),10):
      plt.loglog(dofs[k],errors[k],'-o',lw=1.3, color=palette_list[3], label='_nolegend_')
      plt.loglog(dofs[k][-1],errors[k][-1], marker="o", markersize=10, color=palette_list[3], alpha=alpha, label='_nolegend_')
   plt.loglog(rldofs,rlerrors,marker="^",lw=1.3, color=palette_list[0], label=r'(AM)$^2$R policy')
   plt.loglog(rldofs[-1],rlerrors[-1], marker="^", markersize=10, color=palette_list[0], label='_nolegend_')
   ax4.set_xlabel(r'Degrees of freedom (ndofs$(\mathcal{T}_k)$)', fontsize=22)
   ax4.set_ylabel(r'Global error estimate $(\eta_k)$', fontsize=22)
   ax4.tick_params(axis='x', labelsize=22)
   ax4.tick_params(axis='y', labelsize=22)
   ax4.legend(loc='upper right', prop={'size': 14})
   if save_figs:
      plt.savefig(output_dir+'/'+fig_name_prefix+'_fig4.pdf',format='pdf', bbox_inches='tight')

   ##########
   # make fig5
   ##########
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
   plt.loglog(cumrldofs,rlerrors,marker="^",lw=1.3, color=palette_list[0], label=r'(AM)$^2$R policy')
   # plt.loglog(cumrldofs[-1],rlerrors[-1], marker="o", markersize=10, color=palette_list[0], label='_nolegend_')
   ax5.set_xlabel(r'Cumulative degrees of freedom $(J_k)$', fontsize=22)
   ax5.set_ylabel(r'Global error estimate $(\eta_k)$', fontsize=22)
   ax5.tick_params(axis='x', labelsize=22)
   ax5.tick_params(axis='y', labelsize=22)
   ax5.legend(loc='upper right', prop={'size': 14})
   # for i in range(9,len(errors),10):
   #    print("Cum dofs = ", cumdofs[i][-1])
   # print("cumrldofs = ", cumrldofs[-1])
   if save_figs:
      plt.savefig(output_dir+'/'+fig_name_prefix+'_fig5.pdf',format='pdf', bbox_inches='tight')

   if ex_type == 6 and alpha_fixed: # example 1a MO with fixed alpha; save final error and cumulative dofs
      #dofs = cumrldofs[-1], errors = rlerrors[-1]

      alpha_data_str = alpha_str + ', ' + str(cumrldofs[-1]) + ', ' + str(rlerrors[-1]) + "\n"
      print("writing alpha to file")
      alpha_file = open("const_alpha_data.txt", "a")
      alpha_file.write(alpha_data_str)
      alpha_file.close()

   ##########
   # make fig rlepisode_cost
   ##########
   # ## Make letter-box plot of RL vs expert policies
   plt.figure(figsize=(6,6))
   ax2 = sns.boxenplot(y=costs, width=.6, color=palette_list[3], label='_nolegend_')
   x2 = ax2.get_xlim()
   plt.hlines(rlepisode_cost, x2[0], x2[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy cost')
   y2 = ax2.get_ylim()
   # note: factor 0.1 below makes hatch area have height = 10% of the range of the y2 axis
   plt.fill_between(x2, y2[0] - 0.1 * (y2[1]-y2[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')

   ax2.set_xlabel(r'')
   if minimum_budget_problem:
      ax2.set_ylabel(r'Cost at final step $(\log_2(J_k))$', fontsize=22)
   else:
      ax2.set_ylabel(r'Final error $\left(\log_2(\eta_k)\right)$', fontsize=22)
   ax2.tick_params(axis='y', labelsize=22)
   # lgd = ax2.legend(loc='upper right', prop={'size': 20})
   # letterbox_entry(lgd)
   sns.despine(ax=ax2, bottom=True)
   ax2.tick_params(bottom=False)
   plt.tight_layout()
   if save_figs:
      plt.savefig(output_dir+'/'+fig_name_prefix+'_rlepisode_cost.pdf',format='pdf', bbox_inches='tight')


if ex_type == 4: # example 2c

   ##########
   # make fig rl_vs_two_param
   ##########
   df = pd.read_csv(twopar_file)

   plt.figure(figsize=(6,6))
   ax7 = sns.boxenplot(y=df['costs'], width=.6, color=palette_list[3], label='_nolegend_')
   x7 = ax7.get_xlim()
   plt.hlines(rlepisode_cost, x7[0], x7[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy cost')
   y2 = ax7.get_ylim()
   # plt.fill_between(x7, np.floor(y2[0]), rlepisode_cost, color=palette_list[9], label=r'Apparent performance barrier')
   plt.fill_between(x7, np.floor(y2[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')

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


# plt.show()
