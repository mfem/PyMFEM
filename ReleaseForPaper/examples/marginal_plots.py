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


#########################################
# RL policy df creation or loading
#########################################

file_exists = os.path.exists(output_dir+'/marginals/rl_avg_costs.npy')
if file_exists:
    rl_avg_costs = np.load(output_dir+'/marginals/rl_avg_costs.npy')
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
        print("   cost for angle ", angle, " = ", best_cost)
    avg_cost = sum_of_best_costs/num_angles
    # avg_cost = avg_cost / np.log(2) # convert to log base 2 # don't need to do this - already took log_2
    rl_avg_costs[0] = avg_cost
    print("*** RL avg cost = ", avg_cost, "\n")
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
    num_angles = len(angles)
    num_thetas = len(thetas)
    num_rhos   = len(rhos)
    tpp_avg_costs = np.zeros(num_thetas * num_rhos)
    # if num_angles != 100 or num_thetas != 10 or num_rhos != 10:
    #    print("Something isn't the right size")
    #    exit()
     
        
    i=0
    # temprows = []
    for theta in thetas:
        for rho in rhos:
            sum_of_best_costs = 0.0
            for angle in angles: 
                df_angle = df.loc[(df['theta'] == theta) & (df['rho'] == rho) & (df['angle'] == angle)].set_index('step')
                df_best  = df_angle.iloc[df_angle.index.max()]
                best_cost = df_best['cost']
                sum_of_best_costs += best_cost
            avg_cost = sum_of_best_costs/num_angles
            # avg_cost = avg_cost / np.log(2) # convert to log base 2 # don't need to do this - already took log_2
            tpp_avg_costs[i] = avg_cost
            print(i, " ", theta, " ", rho, " ", avg_cost)
            i += 1
            # temprows.append([theta, rho, avg_cost])
    # trc_df = pd.DataFrame(temprows)
    # trc_df.columns = ['theta', 'rho', 'avgcost']

    np.save(output_dir + '/marginals/tpp_avg_costs.npy', tpp_avg_costs)


#########################################
# expert policy df creation or loading
#########################################

file_exists = os.path.exists(output_dir+'/marginals/exp_avg_costs.npy')
if file_exists:
    exp_avg_costs = np.load(output_dir+'/marginals/exp_avg_costs.npy')
    print("\nLoaded exp_avg_costs.npy\n")
else:
    df = pd.read_csv(output_dir+'/marginals/marginals_exp.csv', index_col=False)   
    angles = df['angle'].unique()
    thetas = df['theta'].unique()
    num_thetas = len(thetas)
    num_angles = len(angles)
    exp_avg_costs = np.zeros(num_thetas)

    # if num_thetas != 100:
    #     print("Should have had 100 theta values... exiting")
    #     exit()
    # if num_angles != 100:
    #     print("Should have had 100 angles... exiting")
    #     exit()

    i=0
    for theta in thetas:
        sum_of_best_costs = 0.0
        for angle in angles: 
            df_angle = df.loc[(df['theta'] == theta) & (df['angle'] == angle)].set_index('step')
            df_best = df_angle.iloc[df_angle.index.max()]
            best_cost = df_best['cost']
            sum_of_best_costs += best_cost
        avg_cost = sum_of_best_costs/num_angles
        # avg_cost = avg_cost / np.log(2) # convert to log base 2 # don't need to do this - already took log_2
        exp_avg_costs[i] = avg_cost
        i += 1
    np.save(output_dir + '/marginals/exp_avg_costs.npy', exp_avg_costs)


case = 101
if case == 101: # marginal distributions for (expert), two policy, and RL polices (100 polices x 100 angles)
 
    f1, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(12,6))  
 
 
 
    ###########################
    # make letterbox plots
    ###########################
 
    num_exp = len(exp_avg_costs)
    num_tpp = len(tpp_avg_costs)
    num_rlp = len(rl_avg_costs)
    num_total = num_exp + num_tpp + num_rlp
    plotdf = pd.DataFrame(np.zeros((num_total,2)))  #  NOTE: assumes 25 RL files
    plotdf.columns = ['policy type', 'avg cost']
    plotdf.iloc[0:num_exp,0] = 'one fixed\n parameter\n "expert" policies'
    plotdf.iloc[0:num_exp,1] = exp_avg_costs
    plotdf.iloc[num_exp:num_exp+num_tpp,0] = 'two fixed\n parameters\n policies'
    plotdf.iloc[num_exp:num_exp+num_tpp,1] = tpp_avg_costs
    plotdf.iloc[num_exp+num_tpp:num_total,0] = 'reinforcement\n learning\n trained policies'
    plotdf.iloc[num_exp+num_tpp:num_total,1] = rl_avg_costs

    print("The best expert    policy was at index ", exp_avg_costs.argmin())
    print("The best two param policy was at index ", tpp_avg_costs.argmin())
 
    ax1 = sns.boxenplot(y='policy type', x='avg cost', data=plotdf, orient="h", palette="Set2")
    # ax1.set_xlabel('Marginal distributions of costs (over angles)')
    ax1.set_title('Marginal distributions on L-shaped type domains',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # ax1.set_xlabel(r'$\mathbb{E}_{\omega}[ \log ( E_{\mathrm{Final}} ) ]$',fontsize=16)
    ax1.set_xlabel(r'$\mathbb{E}_{\omega}[ \log_2 ( E_{\mathrm{Final}} ) ]$',fontsize=16)
    ax1.set_ylabel('')
    # ax1.set_xscale('log')

    plt.savefig(output_dir+'/marginals/fig_101.pdf',format='pdf', bbox_inches='tight')
    filename = output_dir+'/marginals/fig_101.pdf'
    import subprocess
    subprocess.call(["open", filename])
 
    # # Make plot of avg cost as function of theta
    # df = pd.read_csv('../hp_sample_data/most_recent_data/dwf_oct29_all_angles.csv', index_col=False)
    # thetas = df['theta'].unique()
    # x_ax = thetas
    # y_ax = dwf_avg_costs
    # ax2.plot(x_ax, y_ax, 'r', linestyle='-', marker='x', label='dwf avg costs')
    # ax2.set_xlabel('Theta value, expert policy')
    # ax2.set_ylabel('Cost')
 
    #   dwf / dnf: average over angles, this gives an expected cost for each fixed parameter
    # ax3 = sns.boxenplot(y='variable', x = 'Error at final mesh', data=mdf, orient="h", palette="Set2")
 
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



# draw scatter plot of best policies
# plt.scatter(df_bestpolicies['theta'], df_bestpolicies['rho'], c=df_bestpolicies['angle'], cmap='Reds')
# plt.show()
# exit()   


case = 102
if case == 102: # Plot fixed angle action vs. refinement step WITH best two-parameter action overlay

    # numrows = 3
    # numcols = 7

    numrows = 1
    numcols = 5


    plt.rcParams["font.family"] = "Times New Roman"
    plt.rc('text', usetex=True)

    f2, ax2 = plt.subplots(nrows=numrows, ncols=numcols, sharex=True, sharey=True, figsize=(12,3), squeeze=False)  
    # f2, ax2 = plt.subplots(nrows=numrows, ncols=numcols, sharex=True, sharey=True, figsize=(12,6), squeeze=False)  

    df_rl = pd.read_csv(output_dir+'/marginals/marginals_rl.csv', index_col=False)
    angles = df_rl['angle'].unique()

    if numcols == 5:
        angles = np.array([angles[0], angles[5], angles[10], angles[15], angles[20]])

    num_angles = len(angles)
    min_angle = np.round(angles.min(), 2)
    max_angle = np.round(angles.max(), 2)

    df_tpp = pd.read_csv(output_dir+'/marginals/marginals_tpp.csv', index_col=False)   
    assert len(df_tpp['angle'].unique() == num_angles) # check same number of angles tested in both datasets 

    # find two-param policy that achieved best result for each angle
    df_bestpolicies = pd.DataFrame(columns=['angle','theta','rho'])
    for angle in angles:
        df_allperangle  = df_tpp.loc[df_tpp['angle'] == angle][['theta', 'rho', 'cost']]
        df_bestperangle = df_allperangle.iloc[df_allperangle['cost'].argmin()]
        df_bestpolicies = df_bestpolicies.append({'angle': angle, 'theta': df_bestperangle['theta'], 'rho': df_bestperangle['rho']}, ignore_index=True)
    print("Best two parameter policies: \n", df_bestpolicies)


    for i in range(numrows): # num of rows
        for j in range(numcols): 
            index = i*numcols + j
            makelegend=False
            # if index == 7:
            #     makelegend=True
            #     index=6

            rlactions = df_rl[df_rl['angle'] == angles[index]][['theta', 'rho']].to_numpy()
            ax2[i][j].plot(rlactions[:,0],'-o',lw=1.3, c=palette_list[0], label=r'$\theta$ (h parameter)')
            ax2[i][j].plot(rlactions[:,1],'-o',lw=1.3, c=palette_list[1], label=r'$\rho$ (p parameter)')
            # ax2[i][j].plot()
            
            tpp_best_theta = df_bestpolicies[df_bestpolicies['angle'] == angles[index]]['theta'].to_numpy().item() * np.ones(rlactions.shape[0]) 
            tpp_best_rho   = df_bestpolicies[df_bestpolicies['angle'] == angles[index]]['rho'].to_numpy().item()   * np.ones(rlactions.shape[0]) 
            ax2[i][j].plot(tpp_best_theta,'-x',lw=1.3, c=palette_list[0], label=r'$\theta$ (h parameter)')
            ax2[i][j].plot(tpp_best_rho,  '-x',lw=1.3, c=palette_list[1], label=r'$\rho$ (p parameter)')
            ax2[i][j].set_title('angle ' + str(np.round(angles[index]/np.pi,2))+ r'$\pi$')

            if makelegend:
                ax[i][j].legend()


    f2.suptitle('RL vs two-parameter actions for angles ' + str(min_angle) + '-' + str(max_angle) + '\n blue = theta (h); orange = rho (p); outputdir = ' + str(output_dir[-14:]).replace("_"," "), fontsize=18, y=1.2)
    # ttl = f2.suptitle(...)
    # ttl.set_position([.5, 1.01])
    plt.savefig(output_dir+'/marginals/fig_102.pdf',format='pdf', bbox_inches='tight')
    filename = output_dir+'/marginals/fig_102.pdf'
    import subprocess
    subprocess.call(["open", filename])


case = 103
if case == 103: # Alternate comparison of marginals plots, (two param vs RL only)
    plt.figure(figsize=(6,6))
    ax2 = sns.boxenplot(y=tpp_avg_costs, width=.6, color=palette_list[3], label='_nolegend_')
    x2 = ax2.get_xlim()
    plt.hlines(rl_avg_costs, x2[0], x2[1], lw=4, color=palette_list[0], label=r'(AM)$^2$R policy cost')
    y2 = ax2.get_ylim()
    # note: factor 0.1 below makes hatch area have height = 10% of the range of the y2 axis
    # plt.fill_between(x2, y3[0] - 0.1 * (y2[1]-y2[0]), rlepisode_cost, hatch='\\\\\\\\', facecolor=palette_list[9], label=r'Apparent performance barrier')


    ax2.set_ylim(y2[0]-2.0, y2[1]) # forces y axis to leave space for legend at bottom

    ax2.set_xlabel(r'')
    ax2.set_ylabel(r'Average final $\log_2$(error) $\left(\mathbb{E}_{\omega}\left[ \log_2(\eta_k)\right]\right)$', fontsize=22)
    ax2.tick_params(axis='y', labelsize=26)
    lgd = ax2.legend(loc='lower center', prop={'size': 22})
    letterbox_entry(lgd)
    sns.despine(ax=ax2, bottom=True)
    ax2.tick_params(bottom=False)

    plt.tight_layout()

    plt.savefig(output_dir+'/marginals/Ex2c_marginals_letterbox.pdf',format='pdf', bbox_inches='tight')
    filename = output_dir+'/marginals/Ex2c_marginals_letterbox.pdf'
    import subprocess
    subprocess.call(["open", filename])


# rlactions = df_rl[df_rl['angle'] == 0.0][['theta', 'rho']].to_numpy()
# plt.plot(rlactions[:,0],'-o',lw=1.3, label=r'$\theta$ (h parameter)')
# plt.plot(rlactions[:,1],'-o',lw=1.3, label=r'$\rho$ (p parameter)')
# ax2.set_ylabel(r'parameter values in trained (AM)$^2$R policy')
# ax2.legend(loc='upper right')
# ax2.set_xlabel(r'Refinement step')
# plt.savefig(output_dir+'/'+fig_name_prefix+'_fig6.pdf',format='pdf', bbox_inches='tight')


# plt.show()