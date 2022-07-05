"""
    
    EXAMPLE 1: h-refinement policy for L-shaped domain problem

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from prob_envs.Poisson import Poisson
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

"""
    STEP 1: Set parameters
"""

import argparse
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
parser.add_argument('--nbatches', dest='num_batches', type=int, default=150)
parser.add_argument('--demo', dest='demo', action='store_true')
args = parser.parse_args()
print("Parsed options = ", args)
train=args.train
eval=args.eval
save_data=args.savedata
plot_figs=args.plotfigs
save_figs=args.savefigs
nbatches = args.num_batches

restore_policy = False
minimum_budget_problem = True  # mininum budget == error_threshold == minimize dofs

## Configuration for minimum budget problem
prob_config = {
    'mesh_name'         : 'l-shape-benchmark.mesh',
    'problem_type'      : 'lshaped',
    'num_unif_ref'      : 1,
    'order'             : 2,
    'optimization_type' : 'error_threshold', 
    'dof_threshold'     : 5e5,
    'error_threshold'   : 1e-4, # previously 1e-3
    'num_batches'       : nbatches,
    'ee_normalizer'     : 'fixed-order',  # 'fixed-order' for h-ref only or 'variable-order' for hp-ref (or h-ref)
}

## Change to minimum error problem
if not minimum_budget_problem:
    prob_config['optimization_type'] = 'dof_threshold'
    prob_config['dof_threshold']     = 1e4  # previously 5e4

## Neural network configuration
model_config = {
    "fcnet_hiddens"    : [128, 128],
    "fcnet_activation" : "swish",
}

## rllib parameters
config = ppo.DEFAULT_CONFIG.copy()
config['batch_mode'] = 'truncate_episodes'
config['sgd_minibatch_size'] = 100
config['rollout_fragment_length'] = 50
config['num_workers'] = 10
config['train_batch_size'] = config['rollout_fragment_length'] * config['num_workers']
config['num_gpus'] = 0
config['gamma'] = 1.0
config['lr'] = 1e-4
config['seed'] = 4000
config['model'] = model_config

if args.demo:
    config['sgd_minibatch_size'] = 50
    config['rollout_fragment_length'] = 50
    config['num_workers'] = 1
    config['train_batch_size'] = config['rollout_fragment_length'] * config['num_workers']

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
log_dir = homepath + '/capsule/results/logs/'
output_dir_ = homepath + '/capsule/results/output/'

if (restore_policy):
    chkpt_num = nbatches
    # set the path of the checkpoint
    temp_path = 'Example1_2022-04-15_10-55-16'
    checkpoint_dir = log_dir + temp_path
    chkpt_file=checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num)
    output_dir = output_dir_ + temp_path
else:
    timestr = datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    temp_path = 'Example1_' + timestr
    checkpoint_dir = log_dir + temp_path
    output_dir = output_dir_ + temp_path

## Train policy
ray.shutdown()
ray.init(ignore_reinit_error=True,log_to_driver=False)
register_env("my_env", lambda config : Poisson(**prob_config))
trainer = ppo.PPOTrainer(env="my_env", config=config, 
                       logger_creator=custom_log_creator(checkpoint_dir))
env = Poisson(**prob_config)

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
    temp_path = 'Example1_2022-04-15_10-55-16'
    chkpt_num = nbatches
    checkpoint_dir = log_dir + temp_path
    checkpoint_path=checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt > 99 and <1000
    output_dir = output_dir_ + temp_path

if train:
    print_config(checkpoint_dir,prob_config=prob_config, rl_config=rl_config)

"""
    STEP 3: Validation
"""

if eval:

    ## Enact trained policy
    trainer.restore(checkpoint_path)
    rlactions = []
    obs = env.reset()
    done = False
    rlepisode_cost = 0
    rlerrors = [env.global_error]
    rldofs = [env.sum_of_dofs]
    env.trainingmode = False
    while not done:
        action = trainer.compute_single_action(obs,explore=False)
        # action = trainer.compute_single_action(obs,explore=True)
        obs, reward, done, info = env.step(action)
        if not minimum_budget_problem and done:
            break
        rlactions.append(action[0])
        rlepisode_cost -= reward
        print("step = ", env.k)
        print("action = ", action.item())
        print("Num. Elems. = ", env.mesh.GetNE())
        print("episode cost = ", rlepisode_cost)
        rldofs.append(info['num_dofs'])
        rlerrors.append(info['global_error'])
        

    ## Enact AMR policies
    costs = []
    actions = []
    nth = 100
    if args.demo:
        nth = 10
    errors = []
    dofs = []
    for i in range(1, nth):
        print(i)
        action = np.array([i/nth])
        actions.append(action.item())
        print("policy with constant theta = ", action.item())
        env.reset()
        done = False
        episode_cost_tmp = 0
        errors_tmp = [env.global_error]
        dofs_tmp = [env.sum_of_dofs]
        max_steps   = 200
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
        costs.append(episode_cost_tmp)
        errors.append(errors_tmp)
        dofs.append(dofs_tmp)
        print('episode cost = ', episode_cost_tmp)


"""
    STEP 4: Save Data
"""

if (save_data or save_figs):
    mkdir_p(output_dir)
    print_config(output_dir,prob_config=prob_config, rl_config=rl_config)


if save_data:
    import json
    print("Saving data to ", output_dir)
    root_path, _ = os.path.split(checkpoint_path)
    root_path, _ = os.path.split(root_path)
    csv_path = root_path + '/progress.csv'
    df = pd.read_csv(csv_path)
    cost = -df.episode_reward_mean.to_numpy()
    rl_actions = rlactions
    rl_actions.append(rlepisode_cost)
    df1 = pd.DataFrame({'cost':cost})
    df2 = pd.DataFrame({'rlactions':rl_actions,'rldofs':rldofs,'rlerrors':rlerrors})
    df3 = pd.DataFrame({'actions':actions,'costs':costs,'errors':errors,'dofs':dofs})
    filename = output_dir+"/training_data.csv"
    print("Saving training data to: ", filename)    
    df1.to_csv(filename, index=False)
    filename = output_dir+"/rl_data.csv"
    print("Saving RL deployed policy data to: ", filename)    
    df2.to_csv(filename, index=False)
    filename = output_dir+"/expert_data.csv"
    print("Saving deterministic AMR policies data to: ", filename)    
    df3.to_csv(filename, index=False)

"""
    STEP 5: Plots
"""

if plot_figs or save_figs:

    import subprocess
    print("Calling plots.py")
    string_to_call = "python plots.py " + output_dir
    subprocess.call(string_to_call, shell=True)
