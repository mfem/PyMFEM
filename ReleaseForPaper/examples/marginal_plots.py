import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import os.path
import json

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("file_path", type=Path)
# parser.add_argument('--angle_abbrv', type=str, required=False)

args = parser.parse_args()

assert args.file_path.exists() , "The given directory doesn't exist"

output_dir = str(args.file_path)

# if args.angle_abbrv is not None:
#    print("Plotting for angle ", args.angle_abbrv)

if output_dir.find('Example2c') != -1:
   # assert args.angle_abbrv is not None , "Need to provide angle to plots.py for Example 2c"
   print("Loading data from ", output_dir)
   # fig_name_prefix = 'Example2c_angle_' + args.angle_abbrv
   fig_name_prefix = 'Example2c_marginals_'
   ex_type = 4
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

print("*** Check that correct data was loaded here in marginal_plots.py ***")
train_data_file = output_dir+'/training_data.csv'

case = 101
if case == 101: # marginal distributions for (expert), two policy, and RL polices (100 polices x 100 angles)

   f, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(12,6))  



   #########################################
   # RL policy df creation or loading
   #########################################


   file_exists = os.path.exists(output_dir+'/marginals/rl_avg_costs.npy')
   if file_exists:
      tpp_avg_costs = np.load(output_dir+'/marginals/rl_avg_costs.npy')
      print("\nLoaded rl_avg_costs.npy\n")
   else:
      df = pd.read_csv(output_dir+'/marginals/marginals_rl.csv', index_col=False)

      angles = df['angle'].unique()
      num_angles = len(angles)
      # if num_angles != 100:
      #       print("Should have had 100 angles... exiting")
      #       exit()

      rl_avg_costs = np.zeros(1) # 1 indicates only evaluating one RL policy 
      sum_of_best_costs = 0.0
      for angle in angles: 
            df_angle = df.loc[(df['angle'] == angle)].set_index('step')
            df_best  = df_angle.iloc[df_angle.index.max()]
            best_cost = df_best['cost']
            sum_of_best_costs += best_cost
      avg_cost = sum_of_best_costs/num_angles
      avg_cost = avg_cost / np.log(2) # convert to log base 2
      rl_avg_costs[0] = avg_cost
      np.save(output_dir + '/marginals/rl_avg_costs.npy', rl_avg_costs)


   #########################################
   # two param policy df creation or loading #    *** tpp replaced dnf ****
   #########################################

   file_exists = os.path.exists(output_dir+'/marginals/tpp_avg_costs.npy')
   if file_exists:
      tpp_avg_costs = np.load(output_dir+'/marginals/tpp_avg_costs.npy')
      print("\nLoaded tpp_avg_costs.npy\n")
   else:
      df = pd.read_csv(output_dir+'/marginals/marginals_tpp.csv', index_col=False)

      angles = df['angle'].unique()
      thetas = df['theta'].unique()
      rhos   = df['rho'].unique()

      tpp_avg_costs = np.zeros(100)
      num_angles = len(angles)
      num_thetas = len(thetas)
      num_rhos   = len(rhos)
      # if num_angles != 100 or num_thetas != 10 or num_rhos != 10:
      #    print("Something isn't the right size")
      #    exit()
      
      i=0
      for theta in thetas:
         for rho in rhos:
               sum_of_best_costs = 0.0
               for angle in angles: 
                  df_angle = df.loc[(df['theta'] == theta) & (df['rho'] == rho) & (df['angle'] == angle)].set_index('step')
                  df_best  = df_angle.iloc[df_angle.index.max()]
                  best_cost = df_best['cost']
                  sum_of_best_costs += best_cost
               avg_cost = sum_of_best_costs/num_angles
               avg_cost = avg_cost / np.log(2) # convert to log base 2
               tpp_avg_costs[i] = avg_cost
               i += 1
      np.save(output_dir + '/marginals/tpp_avg_costs.npy', tpp_avg_costs)

   ###########################
   # make letterbox plots
   ###########################

   plotdf = pd.DataFrame(np.zeros((225,2)))  #  NOTE: assumes 25 RL files
   plotdf.columns = ['policy type', 'avg cost']
   #  plotdf.iloc[0:100,0] = 'one fixed\n parameter\n "expert" policies'
   #  plotdf.iloc[0:100,1] = dwf_avg_costs
   plotdf.iloc[100:200,0] = 'two fixed\n parameters\n policies'
   plotdf.iloc[100:200,1] = tpp_avg_costs
   plotdf.iloc[200:201,0] = 'reinforcement\n learning\n trained policies'
   plotdf.iloc[200:201,1] = rl_avg_costs

   ax = sns.boxenplot(y='policy type', x='avg cost', data=plotdf, orient="h", palette="Set2")
   # ax.set_xlabel('Marginal distributions of costs (over angles)')
   ax.set_title('Marginal distributions on L-shaped type domains',fontsize=18)
   plt.xticks(fontsize=14)
   plt.yticks(fontsize=14)
   # ax.set_xlabel(r'$\mathbb{E}_{\omega}[ \log ( E_{\mathrm{Final}} ) ]$',fontsize=16)
   ax.set_xlabel(r'$\mathbb{E}_{\omega}[ \log_2 ( E_{\mathrm{Final}} ) ]$',fontsize=16)
   ax.set_ylabel('')
   # ax.set_xscale('log')

   # # Make plot of avg cost as function of theta
   # df = pd.read_csv('../hp_sample_data/most_recent_data/dwf_oct29_all_angles.csv', index_col=False)
   # thetas = df['theta'].unique()
   # x_ax = thetas
   # y_ax = dwf_avg_costs
   # ax.plot(x_ax, y_ax, 'r', linestyle='-', marker='x', label='dwf avg costs')
   # ax.set_xlabel('Theta value, expert policy')
   # ax.set_ylabel('Cost')

   #   dwf / dnf: average over angles, this gives an expected cost for each fixed parameter
   # ax = sns.boxenplot(y='variable', x = 'Error at final mesh', data=mdf, orient="h", palette="Set2")

# rldata_file = output_dir+'/rl_data.csv'
# deterministic_amr_data_file = output_dir+'/deterministic_amr_data.csv'
# if args.angle_abbrv is not None:
#    rldata_file = rldata_file[:-4] + '_angle_' + args.angle_abbrv + '.csv'
#    deterministic_amr_data_file = deterministic_amr_data_file[:-4] + '_angle_' + args.angle_abbrv + '.csv'
#    two_param_data_file = output_dir + '/two_param_amr_data_angle_' + args.angle_abbrv + '.csv'


# df = pd.read_csv(train_data_file, index_col=0)
# cost = df.index.to_numpy()

# if ex_type == 4: # example 2c
#    df = pd.read_csv(rldata_file)
#    rlepisode_cost = df['rlepisode_cost'][0]
#    rlactions = df[{'theta', 'rho'}].to_numpy()
#    rlactions = rlactions[:-1,:] # drop last row, which was padded
#    rldofs =  df['rldofs'].to_numpy()
#    rlerrors =  df['rlerrors'].to_numpy()
# else:
#    df = pd.read_csv(rldata_file, index_col=0)
#    rlactions = df.index.to_numpy()
#    rlepisode_cost = rlactions[-1]
#    rlactions = rlactions[:-1]
#    rldofs = df.iloc[:, 0].to_numpy()
#    rlerrors = df.iloc[:, 1].to_numpy()


# df = pd.read_csv(deterministic_amr_data_file, index_col=0)

# actions = df.index.to_numpy()
# costs = df.iloc[:, 0].to_numpy()
# errors = df.iloc[:, 1].to_numpy()
# dofs = df.iloc[:, 2].to_numpy()

# for i in range(len(dofs)):
#    dofs[i] = np.array(eval(dofs[i]))
#    errors[i] = np.array(eval(errors[i]))

plt.show()