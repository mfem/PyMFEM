"""
    
    EXAMPLE 2-c multi-objective: h-refinement policy for L-shaped domain with missing angle problem with multi-objective cost function

"""

import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from prob_envs.MultiObjectivePoisson import hp_Angle_MultiObjPoisson
import numpy as np
import time
import seaborn as sns

from ray.tune.logger import UnifiedLogger
from datetime import datetime
import json 

def print_config(dir, prob_config = None, rl_config = None):
    if (prob_config is not None):
        with open(dir+"/prob_config.json", 'w') as f: 
            json.dump(prob_config,f)
            # for key, value in prob_config.items(): 
                # f.write('%s:%s\n' % (key, value))

    if (rl_config is not None):
        with open(dir+"/rl_config.json", 'w') as f: 
            json.dump(rl_config,f)

def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise


def custom_log_creator(custom_path):

    logdir_prefix = custom_path
    def logger_creator(config):

        if not os.path.exists(custom_path):
            os.makedirs(custom_path)
        logdir = logdir_prefix
        return UnifiedLogger(config, logdir, loggers=None)

    return logger_creator


sns.set()
# sns.set_context("notebook")
# sns.set_context("paper")
# sns.set_context("paper", font_scale=1.5)
sns.set_context("talk", font_scale=3)
# sns.set_style("ticks")
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
# sns.set_theme(style="white", rc=custom_params)

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

"""
    STEP 1: Set parameters
"""

parser = argparse.ArgumentParser()
parser.add_argument('--train', default=True, action='store_true')
parser.add_argument('--no-train', dest='train', action='store_false')
parser.add_argument('--eval', default=True, action='store_true')
parser.add_argument('--no-eval', dest='eval', action='store_false')
parser.add_argument('--savedata', default=True, action='store_true')
parser.add_argument('--no-savedata', dest='savedata', action='store_false')
parser.add_argument('--plotfigs', default=True, action='store_true')
parser.add_argument('--no-plotfigs', dest='plotfigs', action='store_false')
parser.add_argument('--savefigs', default=True, action='store_true')
parser.add_argument('--no-savefigs', dest='savefigs', action='store_false')
parser.add_argument('--savemesh', default=False, action='store_true')

parser.add_argument('--angle', default = 0.5*np.pi, type = float)

parser.add_argument('--observe_alpha', default = True, action='store_true')
parser.add_argument('--no_observe_alpha', dest='observe_alpha', action='store_false')
parser.add_argument('--alpha', default = 0.5, type = float) # alpha for cost function if no_observe_alpha, also the alpha used for evaluation

parser.add_argument('--observe_budget', default = True, action='store_true')
parser.add_argument('--no_observe_budget', dest='observe_budget', action='store_false')


args = parser.parse_args()
print("Parsed options = ", args)
train=args.train
eval=args.eval
save_data=args.savedata
plot_figs=args.plotfigs
save_figs=args.savefigs

restore_policy = False
nbatches = 300

## Configuration for multi objective problem
prob_config = {
    'mesh_name'         : 'l-shape-benchmark.mesh',
    'problem_type'      : 'lshaped',
    'num_unif_ref'      : 1,
    'order'             : 2,
    'optimization_type' : 'multi_objective', 
    'dof_threshold'     : 1e5,
    'alpha'             : args.alpha,
    'observe_alpha'     : args.observe_alpha,
    'observe_budget'    : args.observe_budget,
    'num_iterations'    : 15,
    'num_batches'       : nbatches
}

mesh_abbrv = prob_config['problem_type']

## Change to minimum error or minimum dof problem
if prob_config['optimization_type'] == 'error_threshold': # minimum dof
    prob_config['error_threshold']   = 1e-4

elif prob_config['optimization_type'] == 'dof_threshold': #minimum error
    prob_config['dof_threshold']     = 1e4
    prob_config['error_threshold']   = 1e-4

## Neural network configuration
model_config = {
    "fcnet_hiddens"    : [128, 128],
    "fcnet_activation" : "swish",
}

## rllib parameters
config = ppo.DEFAULT_CONFIG.copy()
config['batch_mode'] = 'truncate_episodes'
# config['batch_mode'] = 'complete_episodes'
config['sgd_minibatch_size'] = 100
config['rollout_fragment_length'] = 50
config['num_workers'] = 10
config['train_batch_size'] = config['rollout_fragment_length'] * config['num_workers']
config['num_gpus'] = 0
config['gamma'] = 1.0
config['lr'] = 5e-6
config['seed'] = 4000
config['model'] = model_config

# for limited printing of rl_config
rl_config = {
    'batch_mode'                : config['batch_mode'],
    'sgd_minibatch_size'        : config['sgd_minibatch_size'],
    'rollout_fragment_length'   : config['rollout_fragment_length'],
    'num_workers'               : config['num_workers'],
    'train_batch_size'          : config['train_batch_size'],
    'num_gpus'                  : config['num_gpus'],
    'gamma'                     : config['gamma'],
    'lr'                        : config['lr'],
    'seed'                      : config['seed'],
    'model'                     : config['model'],
}

"""
    STEP 2: Training
"""

homepath = os.path.expanduser("~")
log_dir = os.getcwd() + '/logs/'
output_dir_ = os.getcwd() + '/output/'

if (restore_policy):
    chkpt_num = nbatches
    # set the path of the checkpoint
    temp_path = 'Example1a_2022-04-15_10-55-16'
    checkpoint_dir = log_dir + temp_path
    chkpt_file=checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num)
    output_dir = output_dir_ + temp_path
else:
    timestr = datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    angle_abbrv = "{:.2f}".format(np.round(args.angle,2)) # to keep filename length short

    if prob_config['optimization_type'] == 'multi_objective':
        if args.observe_alpha == True:
            temp_path = 'Example2c_MO_ang_' + angle_abbrv + "_" + timestr
        else: 
            alpha_str = str(args.alpha).replace('.','_') + '_'
            temp_path = 'Example2c_MO_alpha' + alpha_str + "_ang_" + angle_abbrv + "_" + timestr
        
    else:
        temp_path = 'Example2c_' + "_ang_" + angle_abbrv + "_" + timestr

    checkpoint_dir = log_dir + temp_path
    output_dir = output_dir_ + temp_path   

## Train policy
ray.shutdown()
ray.init(ignore_reinit_error=True)
register_env("my_env", lambda config : hp_Angle_MultiObjPoisson(**prob_config))
trainer = ppo.PPOTrainer(env="my_env", config=config, 
                       logger_creator=custom_log_creator(checkpoint_dir))
env = hp_Angle_MultiObjPoisson(**prob_config)

if (restore_policy):
    trainer.restore(chkpt_file)

if train:
    env.trainingmode = True
    for n in range(nbatches):
        print("training batch %d of %d batches" % (n+1,nbatches))
        result = trainer.train()
        episode_score = -result["episode_reward_mean"]
        episode_length = result["episode_len_mean"]
        print ("Mean episode cost:   %.3f" % episode_score)
        print ("Mean episode length: %.3f" % episode_length)
        checkpoint_path = trainer.save(checkpoint_dir)
        print(checkpoint_path)

if eval and not train:
    if prob_config['optimization_type'] == 'multi_objective':
#        temp_path = 'Example1a_MO_2022-07-11_06-56-00'
        temp_path = 'Example1a_MO_2022-07-13_11-56-00'
    else:
        temp_path = 'Example1a_2022-04-15_10-55-16'

    chkpt_num = nbatches
    checkpoint_dir = log_dir + temp_path
    if chkpt_num < 100:
        checkpoint_path=checkpoint_dir+'/checkpoint_0000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt < 100
    elif chkpt_num > 100:
        checkpoint_path = checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt > 99 and <1000
    else:
        print("error, cannot load policy to evaluate")
    output_dir = output_dir_ + temp_path

if train:
    print_config(checkpoint_dir, prob_config=prob_config, rl_config=rl_config)

"""
    STEP 3: Validation
"""

if eval:

    if prob_config['problem_type'] == 'lshaped':
        env.set_angle(args.angle) # note: set_angle(...) redefines self.initial_mesh
        print("*** Set angle for eval to  ", args.angle)

    ## Enact trained policy
    trainer.restore(checkpoint_path)
    cols = ['theta', 'rho']
    rlactions = pd.DataFrame(columns=cols)
    obs = env.reset(random_angle=False, new_alpha = False)
    done = False
    rlepisode_cost = 0
    rlerrors = [env.global_error]
    rldofs = [env.sum_of_dofs]
    env.trainingmode = False

    if args.savemesh:
        mkdir_p(output_dir+"/meshes_and_gfs/")
        env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + mesh_abbrv + "_angle_" + str(angle_abbrv) + '_initial.mesh')
        print("==> Saved initial mesh to ", output_dir, "/meshes_and_gfs/")

    while not done:
        action = trainer.compute_single_action(obs, explore=False)
        obs, reward, done, info = env.step(action)
        if prob_config['optimization_type'] == 'dof_threshold' and done:
            break
        rlactions = rlactions.append({'theta': action[0], 'rho':action[1]}, ignore_index=True)
        rlepisode_cost -= reward

        # print episode results
        print("step = ", env.k)
        print("theta action = ", action[0].item())
        print("rho   action = ", action[1].item())
        print("Num. Elems. = ", env.mesh.GetNE())
        print("episode cost = ", rlepisode_cost)

        rldofs.append(info['num_dofs'])
        rlerrors.append(info['global_error'])

        if args.savemesh:
            mkdir_p(output_dir+"/meshes_and_gfs/")
            gfname = output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + mesh_abbrv + "_angle_" + str(angle_abbrv) + '.gf'
            env.RenderHPmesh(gfname=gfname)
            env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + mesh_abbrv + "_angle_" + str(angle_abbrv) + '.mesh')

    if train == False and prob_config['optimization_type'] == 'multi_objective' and args.observe_alpha == True:
        # save final errors in file
        file_name = "alpha_policy_data_2c.txt"
        file_location = output_dir_ + file_name
        file = open(file_location, "a")

        cum_rldofs = np.cumsum(rldofs)
        file_string = str(args.alpha) + ", " + str(args.angle) + ", " + str(cum_rldofs[-1]) + ", " + str(rlerrors[-1]) + "\n"
        file.write(file_string)
        file.close()
        

    #############
    # Enact fixed "two parameter" AMR policies
    ########### 
    tp_costs = []
    tp_nth = 10
    # tp_actions = np.zeros(((tp_nth-1)**2,2)) # exclude 0.0, 1.0 as actions
    # tp_actions = np.zeros(((tp_nth+1)**2,2)) # include 0.0, 1.0 as actions
    tp_actions = np.zeros(((tp_nth)**2,2)) # include 0.0 but not 1.0 as actions
    # tp_actions = np.zeros((1,2)) # include only theta = [one fixed value], rho = [one fixed value]
    tp_errors = []
    tp_dofs = []
    index_count = 0

    # for theta in range(1, tp_nth):      # exclude 0.0, 1.0 as actions
    #     for rho in range(1, tp_nth):    # exclude 0.0, 1.0 as actions
    # for theta in range(0, tp_nth+1):    # include 0.0, 1.0 as actions
    #     for rho in range(0, tp_nth+1):  # include 0.0, 1.0 as actions
    for theta in range(0, tp_nth):        # include 0.0 but not 1.0 as actions
        for rho in range(0, tp_nth):      # include 0.0 but not 1.0 as actions
    # for theta in range(6, 7):        # include only theta = 0.6, rho = 0.5
    #     for rho in range(5, 6):      # include only theta = 0.6, rho = 0.5
    # for theta in range(6, 7):        # include only theta = 0.6, rho = 0.4
    #     for rho in range(4, 5):      # include only theta = 0.6, rho = 0.4
    # for theta in range(6, 7):        # include only theta = 0.6, rho = 0.3
    #     for rho in range(3, 4):      # include only theta = 0.6, rho = 0.3
            tp_actions[index_count] = np.array([theta/tp_nth, rho/tp_nth]) # note 1D action space
            if theta/tp_nth == 1 and rho/tp_nth == 1: # avoid some linear algerbra error if action is [1,1]
                tp_actions[index_count] = np.array([0.99, 0.99])
            action = tp_actions[index_count]
            print("two param action = ", tp_actions[index_count])
            index_count += 1
            obs = env.reset(random_angle=False, new_alpha = False)
            done = False
            episode_cost_tmp = 0
            errors_tmp = [env.global_error]
            dofs_tmp = [env.sum_of_dofs]
            max_steps   = 40 # or:  num_steps_of_RL_policy
            steps_taken = 0
            while not done:
                _, reward, done, info = env.step(action)
                if not minimum_budget_problem and done:
                    break
                if steps_taken > max_steps:
                    print("*** BREAKING EARLY - fixed action exceeded max step threshold of ", max_steps, "steps.")
                    break
                else:
                    steps_taken += 1
                episode_cost_tmp -= reward
                errors_tmp.append(info['global_error'])
                dofs_tmp.append(info['num_dofs'])
            tp_costs.append(episode_cost_tmp)
            tp_errors.append(errors_tmp)
            tp_dofs.append(dofs_tmp)
            print('two param episode cost = ', episode_cost_tmp)
            # if theta == 5 and rho == 5:
            #     env.render()
            #     env.RenderHPmesh()
            #     print("\nRendering two parmaeter policy ", tp_actions[index_count-1], "\n")
            if save_mesh and prob_config['mesh_name'] == 'fichera.mesh':
                mkdir_p(output_dir+"/meshes_and_gfs/")
                gfname = output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + mesh_abbrv + "_angle_" + str(angle_abbrv) + '_tpp.gf'
                env.RenderHPmesh(gfname=gfname)
                env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + mesh_abbrv + "_angle_" + str(angle_abbrv) + '_tpp.mesh')


"""
    STEP 4: Save Data
"""

if (save_data or save_figs):
    mkdir_p(output_dir)
    print_config(output_dir, prob_config=prob_config, rl_config=rl_config)


if save_data:
    print("Saving data to ", output_dir)
    root_path, _ = os.path.split(checkpoint_path)
    root_path, _ = os.path.split(root_path)
    csv_path = root_path + '/progress.csv'
    df = pd.read_csv(csv_path)
    cost = -df.episode_reward_mean.to_numpy()

    #### training data df
    df1 = pd.DataFrame({'cost':cost})
    filename = output_dir+"/training_data.csv"
    print("Saving training data to: ", filename)    
    df1.to_csv(filename, index=False)

    #### rl action df
    # pad rlactions with one extra row to enable df2 creation
    rlactions = rlactions.append({'theta': 999, 'rho': 999}, ignore_index=True)
    df2 = pd.DataFrame({'theta':rlactions.iloc[:,0], 'rho':rlactions.iloc[:,1], 'rldofs':rldofs,'rlerrors':rlerrors})
    # save a single value in every row of new column 'rlepisode_cost'
    df2['rlepisode_cost'] = rlepisode_cost 
    filename = output_dir+"/rl_data_" + mesh_abbrv + "_angle_" + angle_abbrv + ".csv"
    print("\nSaving RL deployed policy data to: ", filename)  
    df2.to_csv(filename, index=False)

    #### expert policy df
    
    if prob_config['problem_type'] == 'lshaped' and False:
        df3 = pd.DataFrame({'actions':actions,'costs':costs,'errors':errors,'dofs':dofs})
        filename = output_dir+"/expert_data_" + mesh_abbrv + "_angle_" + angle_abbrv + ".csv"
        print("\nSaving expert AMR policies data to: ", filename)    
        df3.to_csv(filename, index=False)

    #### two param policy df
    df4 = pd.DataFrame({'theta':tp_actions[:,0], 'rho':tp_actions[:,1],'costs':tp_costs,'errors':tp_errors,'dofs':tp_dofs})
    filename = output_dir + "/tpp_data_" + mesh_abbrv + "_angle_" + angle_abbrv + ".csv"
    print("\nSaving two parameter AMR policies data to: ", filename)    
    df4.to_csv(filename, index=False)

"""
    STEP 5: Plots
"""

if plot_figs or save_figs:

    import subprocess
    print("Calling plots.py")
    string_to_call = "python plots.py " + output_dir + " --angle " + angle_abbrv + " --mesh " + mesh_abbrv
    subprocess.call(string_to_call, shell=True)

    # print name of output_dir to file for plotting with slurm scritps
    file = open("output_dir4plots.txt","a")
    file.write("\n" + output_dir)
    file.close() 
