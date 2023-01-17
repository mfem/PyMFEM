"""
    
    EXAMPLE 1-A multi-objective: h-refinement policy for L-shaped domain problem with multi-objective cost function

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from prob_envs.MultiObjectivePoisson import MultiObjPoisson
from prob_envs.MultiObjectivePoisson import ADF_MultiObjPoisson
from MO_eval import *
import numpy as np
import time
import seaborn as sns
from gym import spaces

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

# define global variables for tau and delta
@ray.remote
class ADF_Parameters:
    #def __init__(self):
      

    def set_initial_parameters(self, tau_init, tau_inc, delta_warm, delta_anneal):
        self.tau = tau_init
        self.tau_max = tau_init
        self.tau_init = tau_init
        self.tau_inc = tau_inc

        self.delta = delta_warm
        self.delta_anneal = delta_anneal

        # keep track of largest tau used so far
        self.last_tau = self.tau
    def get_delta(self):
        return self.delta
    def update_delta(self):
        self.delta = self.delta #self.delta_anneal

    def update_tau(self):
        self.tau = self.last_tau - self.tau_inc
        self.last_tau = self.tau
    def get_tau(self):
        return self.tau
    def get_tau_init(self):
        return self.tau_init
    def set_tau(self, new_tau):
        self.tau = new_tau
    def set_rand_tau(self):
        # set tau back to a tau we have visited (or almost visited) before
        # this random tau is uniformly distributed between the most recently trained tau
        # and the initial tau, tau_max.
        self.tau = np.random.uniform(low = self.last_tau, high = self.tau_max)

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
parser.add_argument('--savemesh', default=False, action='store_true')

parser.add_argument('--observe_alpha', default = True, action='store_true')
parser.add_argument('--no_observe_alpha', dest='observe_alpha', action='store_false')
parser.add_argument('--alpha', default = 0.5, type = float) # alpha for cost function if no_observe_alpha, also the alpha used for evaluation
parser.add_argument('--eval_alphas', default = False, action='store_true') # evaluate trained policy on alphas in {0.1, 0.2, ..., 0.9}

parser.add_argument('--observe_budget', default = True, action='store_true')
parser.add_argument('--no_observe_budget', dest='observe_budget', action='store_false')


args = parser.parse_args()
print("Parsed options = ", args)

train       = args.train
eval        = args.eval
save_data   = args.savedata
plot_figs   = args.plotfigs
save_figs   = args.savefigs

restore_policy = False
nbatches       = 50

## Configuration for multi objective problem
prob_config = {
    'mesh_name'         : 'l-shape-benchmark.mesh', #'star.mesh',
    'problem_type'      : 'lshaped', #'default',
    'num_unif_ref'      : 1,
    'order'             : 2,
    'optimization_type' : 'multi_objective', 
    'cost_function'     : 'ADF', # ADF for annealing desirability function, alpha for linear combination of costs
    'dof_threshold'     : 5e4,
    'alpha'             : args.alpha,
    'observe_alpha'     : args.observe_alpha,
    'observe_budget'    : args.observe_budget,
    'num_iterations'    : 30,
    'num_batches'       : nbatches,
    'tau_min'           : np.log2(1e-6), #np.log2(1e-2), # np.log2(1e-3), #np.log2(1e-4),
    'M_warm'            : 50, # number of batches in warming phase
    'M_anneal'          : 30,  #5, # number of batches per tau in annealing phase
    'N_anneal'          : 3,  #20,  # number of target errors (tau) to train on (not counting intitial tau)
    'M_retrain'         : 0,  #2,  # number of batches for retraining before each new tau
    'batch_size'        : 100 # number of episodes per batch
}

# if using ADF algorithm, the number of batches is defined by the values
# of M_warm, M_anneal, and N_anneal
if prob_config['cost_function'] == 'ADF':
    M_warm   = prob_config['M_warm']; 
    M_anneal = prob_config['M_anneal']; 
    N_anneal = prob_config['N_anneal'];
    M_retrain = prob_config['M_retrain']
    nbatches = M_warm + (M_anneal + M_retrain) * N_anneal # the +M_retrain corresponds to an extra training batch to retrain on taus we've already seen
    prob_config['num_batches'] = nbatches

    ADF_bool   = True 
    alpha_bool = False 
elif prob_config['cost_function'] == 'alpha':
    ADF_bool   = False 
    alpha_bool = True

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
# config['batch_mode']            = 'truncate_episodes'
config['batch_mode']              = 'complete_episodes'
config['sgd_minibatch_size']      = 100
config['rollout_fragment_length'] = 50
config['num_workers']             = 10
config['train_batch_size']        = config['rollout_fragment_length'] * config['num_workers']
config['num_gpus']                = 0
config['gamma']                   = 1.0
config['lr']                      = 1e-5
config['seed']                    = 4000
config['model']                   = model_config
config['entropy_coeff']           = 0.05

# set up learning rate schedule
# lr_schedule = ray.rllib.utils.schedules.linear_schedule.LinearSchedule(initial_p = 10**-4, final_p = 10**-5, schedule_timesteps = nbatches, framework = 'tf') 
#config['lr_schedule'] = lr_schedule

# # set up epsilon greedy exploration
schedule_object = ray.rllib.utils.schedules.constant_schedule.ConstantSchedule(value = 0.04, framework = 'tf')
action_space = spaces.Box(low=0.0, high=0.999, shape=(1,), dtype=np.float32)
config['explore'] = True
# config['exploration_config'] = {
#    # Exploration sub-class by name or full path to module+class
#    # (e.g. “ray.rllib.utils.exploration.epsilon_greedy.EpsilonGreedy”)
#    'type'               : 'EpsilonGreedy',
#    # Parameters for the Exploration class' constructor:
#    'epsilon_schedule'   : schedule_object,
#    'action_space'       : action_space
#    }

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

    if prob_config['optimization_type'] == 'multi_objective':
        if ADF_bool:
            temp_path = 'Example1a_ADF_' + timestr

        elif alpha_bool:
            if args.observe_alpha == True:
                temp_path = 'Example1a_MO_' + timestr
            else: 
                alpha_str = str(args.alpha).replace('.','_') 
                temp_path = 'Example1a_MO_alpha' + alpha_str + '_' + timestr
        
    else:
        temp_path = 'Example1a_' + timestr

    checkpoint_dir = log_dir + temp_path
    output_dir = output_dir_ + temp_path   

## Train policy
ray.shutdown()
ray.init(ignore_reinit_error=True)

if alpha_bool:
    register_env("my_env", lambda config : MultiObjPoisson(**prob_config))
    trainer = ppo.PPOTrainer(env="my_env", config=config, 
                           logger_creator=custom_log_creator(checkpoint_dir))
    env = MultiObjPoisson(**prob_config)

elif ADF_bool:
    # initialize parameter class
    ADF_params = ADF_Parameters.options(name = "parameters").remote()

    register_env("my_env", lambda config : ADF_MultiObjPoisson(**prob_config))
    trainer = ppo.PPOTrainer(env="my_env", config=config, 
                           logger_creator=custom_log_creator(checkpoint_dir))
    env = ADF_MultiObjPoisson(**prob_config)
    print("env = ", env)
    print("prob config = ", prob_config)

else:
    print("Error: invalid choice of cost function")


if (restore_policy):
   trainer.restore(chkpt_file)

# initialize tau and delta in ADF_Parameter class
if ADF_bool:
    # tau_max, tau_step, delta_warm, delta_anneal = env.FindParameters()
    # print("found params: ", tau_max, tau_step, delta_warm, delta_anneal)
    # 

    # print("WARNING: setting initial tau as tau_min in line below; should be tau_max")
    # ADF_params.set_initial_parameters.remote(prob_config['tau_min'], tau_step, delta_warm, delta_anneal)

    print("\n\n\n ****** WARNING: hard coding initial parameters \n\n") # ***************

    tau_max = np.log2(5e-2)
    tau_step = np.log2(5e-2) - np.log2(1e-2)
    delta_warm = 1
    delta_anneal = 1

    print("assiging params: ", tau_max, tau_step, delta_warm, delta_anneal)
    ADF_params.set_initial_parameters.remote(tau_max, tau_step, delta_warm, delta_anneal)
    
    # print("assiging params: ", np.log2(1e-2), (np.log2(1e-2)-np.log2(1e-3)), delta_warm, delta_anneal)
    # ADF_params.set_initial_parameters.remote(np.log2(1e-2), (np.log2(1e-2)-np.log2(1e-3)), delta_warm, delta_anneal)
    # print("assiging params: ", np.log2(1e-3), (np.log2(1e-3)-np.log2(1e-4)), delta_warm, delta_anneal)
    # ADF_params.set_initial_parameters.remote(np.log2(1e-3), (np.log2(1e-3)-np.log2(1e-4)), delta_warm, delta_anneal)

if train:
    env.trainingmode = True
    MO_eval_loss = []

    for n in range(nbatches):
        # print("  ==> commented out code below to update tau during training <==")
        if ADF_bool:
            # uncomment for ADF approach 
            
            if n == M_warm:
                # we've finished the warming phase; update tau and delta
                ADF_params.update_tau.remote()
                ADF_params.update_delta.remote()

                print("\n\nUpdated delta to {}".format(ray.get(ADF_params.get_delta.remote())))

            elif n> M_warm and (n - M_warm) % (M_anneal + M_retrain) == 0:
                ADF_params.update_tau.remote()

                print("\n\nUpdated tau to {}".format(ray.get(ADF_params.get_tau.remote())))
            elif n > M_warm and (n - M_warm + 1) % (M_anneal + M_retrain) == 0:
                # retrain on old tau
                ADF_params.set_rand_tau.remote()
                print("\n\nretraining on randomly sampled tau = {}".format(ray.get(ADF_params.get_tau.remote())))
            
        print("training batch %d of %d batches" % (n+1, nbatches))
        print("tau = {}".format(ray.get(ADF_params.get_tau.remote())))
        result = trainer.train()
        episode_score = -result["episode_reward_mean"]
        episode_length = result["episode_len_mean"]
        print ("Mean episode cost:   %.3f" % episode_score)

        checkpoint_path = trainer.save(checkpoint_dir)
        print(checkpoint_path)

if eval and not train:
    MO_eval_loss = []
    if prob_config['optimization_type'] == 'multi_objective':
        # temp_path = 'Example1a_MO_2022-07-11_06-56-00'    # experiment set 8
        # temp_path = 'Example1a_MO_2022-07-28_08-10-58'    # experiment set 17
        # temp_path = 'Example1a_MO_2022-07-28_12-46-24'    # experiment set 18
        # temp_path = 'Example1a_MO_2022-07-28_14-22-54'    # experiment set 19
        # temp_path = 'Example1a_MO_2022-07-29_07-23-24'    # experiment set 20
        # temp_path = 'Example1a_MO_2022-07-29_13-16-49'    # experiment set 22
        # temp_path = 'Example1a_MO_2022-07-29_08-54-14'    # experiment set 21
        # temp_path = 'Example1a_MO_2022-08-01_07-46-52'    # experiment set 23
        # temp_path = 'Example1a_MO_2022-08-01_09-30-28'    # experiment set 24
        temp_path = 'Example1a_ADF_2023-01-16_00-04-27'
    else:
        print("*** need to set path for eval *** exiting")
        exit()
        # temp_path = 'Example1a_2022-04-15_10-55-16'

    chkpt_num = nbatches
    checkpoint_dir = log_dir + temp_path
    if chkpt_num < 100:
        checkpoint_path=checkpoint_dir+'/checkpoint_0000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt < 100
    elif chkpt_num >= 100 and chkpt_num <1000:
        checkpoint_path = checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt > 99 and <1000
    elif chkpt_num >=1000:
        checkpoint_path = checkpoint_dir+'/checkpoint_00'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt > 999 and <10000
    else:
        print("error, cannot load policy to evaluate")
    output_dir = output_dir_ + temp_path

if train:
    print_config(checkpoint_dir, prob_config=prob_config, rl_config=rl_config)

print("****\n WARNING: see above in code where initial tau is set as tau_min; should be tau_max;\n  also need to adjust in FindParameters \n****")

"""
    STEP 3: Validation
"""

print("\n\n\n Validation phase - RL policy tests \n\n\n")

# determine which, if any, alphas to evaluate
if alpha_bool:
    if eval and not args.eval_alphas:
        params_to_eval = [args.alpha]

    elif eval and args.eval_alphas:
        params_to_eval = 0.1*np.array(range(1, 10, 1))
    else:
        params_to_eval = []

# determine which errors/taus to evaluate on (tau is the 2log of the target error)
elif ADF_bool:
    if eval:
        tau_max = np.log2(5e-2)
        tau_step = np.log2(5e-2) - np.log2(1e-2)
        params_to_eval = np.array([tau_max - i * tau_step for i in range(5)])
        # params_to_eval = np.array([prob_config['tau_min']])
        print("\n *** params to eval =", params_to_eval, "\n")
        # min_error = env.tau_min
        # max_error = tau_max
        # step_size = (max_error - min_error)/10
        # params_to_eval = min_error + step_size*np.array(range(0,11,1))
    else:
        params_to_eval = []

# create file location for saving evaluation results
mkdir_p(output_dir)
if alpha_bool:
    file_name = "/alpha_policy_data_1a.txt"
elif ADF_bool:
    file_name = "/ADF_policy_data_1a.txt"
file_location = output_dir + file_name

if os.path.exists(file_location):
    os.remove(file_location)
    print("==> Removed old ADF policy data file")

# boolean to keep track of saving mesh, so we don't waste time doing this more than once
if args.savemesh:
    saved_initial_mesh = False;
    mkdir_p(output_dir+"/meshes_and_gfs/")
    env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + prob_config['problem_type'] + '_initial.mesh')

    print("==> Saved initial mesh to ", output_dir, "/meshes_and_gfs/")
    saved_initial_mesh = True;

for param in params_to_eval:
    ## Enact trained policy
    trainer.restore(checkpoint_path) 
    rlactions = []

    if alpha_bool:
        # set alpha
        env.alpha = param
        alpha_str = str(param).replace('.','_') 
        obs = env.reset(new_alpha = False)
    elif ADF_bool:
        # set tau
        ADF_params.set_tau.remote(param)
        tau_str = str(param).replace('.','_')
        obs = env.reset()
    
    done = False
    rlepisode_cost = 0
    rlerrors = [env.global_error]
    rldofs = [env.sum_of_dofs]
    env.trainingmode = False

    while not done:
        # action = trainer.compute_single_action(obs, explore = False)
        action = trainer.compute_single_action(obs,explore=True)
        obs, reward, done, info = env.step(action)
        if prob_config['optimization_type'] == 'dof_threshold' and done:
            break
        rlactions.append(action[0])
        rlepisode_cost -= reward
        print("RL step =", env.k, 
              "axn =", np.round(float(action.item()),3),
              "dof=", info['num_dofs'], 
              # dofs at this step = env.mesh.GetNE(),
              "err=", np.round(info['global_error'],6),
              "cum cost=", rlepisode_cost)
        rldofs.append(info['num_dofs']) 
        rlerrors.append(info['global_error'])

    print("")

    if args.savemesh:
        #gfname = output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + prob_config['problem_type'] + "_alpha_" + alpha_str + '.gf'
        env.render()
        if alpha_bool:
            env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + prob_config['problem_type'] + "_alpha_" + alpha_str + '.mesh')
        elif ADF_bool:
            env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + prob_config['problem_type'] + "_tau_" + tau_str + '.mesh')

    # save final errors in file for each alpha/tau
    cum_rldofs = np.cumsum(rldofs)
    if alpha_bool:
        file_string = str(env.alpha) + ", " + str(cum_rldofs[-1]) + ", " + str(rlerrors[-1]) + "\n"
    elif ADF_bool:
        file_string = str(ray.get(ADF_params.get_tau.remote())) + ", " + str(cum_rldofs[-1]) + ", " + str(rlerrors[-1]) + ", " + str(env.k) + "\n"
    file = open(file_location, "a")
    file.write(file_string)
    file.close()
            
print("\n\n\n Validation phase - Fixed theta policy tests \n\n\n")

## Enact AMR policies
costs = []
actions = []
nth = 100 # number of fixed policies to try
errors = []
dofs = []
for i in range(1, nth):
    print(i)
    action = np.array([i/nth])
    actions.append(action.item())
    print("action = ", action.item())
    if alpha_bool:
        env.reset(new_alpha = False)
    elif ADF_bool:
        env.reset()
    done = False
    episode_cost_tmp = 0
    errors_tmp = [env.global_error]
    dofs_tmp = [env.sum_of_dofs]
    # max_steps   = prob_config['num_iterations']
    steps_taken = 0
    
    while not done:
        _, reward, done, info = env.step(action)
        if prob_config['optimization_type'] == 'dof_threshold' and done:
            break
        if steps_taken > prob_config['num_iterations']:
            print("*** fixed action exceeded max step threshold of ", prob_config['num_iterations'], "steps.")
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
    print_config(output_dir, prob_config=prob_config, rl_config=rl_config)


if save_data:
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
    df4 = pd.DataFrame({'MO_eval_loss':MO_eval_loss})
    filename = output_dir+"/training_data.csv"
    print("Saving training data to: ", filename)    
    df1.to_csv(filename, index=False)
    filename = output_dir+"/rl_data.csv"
    print("Saving RL deployed policy data to: ", filename)    
    df2.to_csv(filename, index=False)
    filename = output_dir+"/expert_data.csv"
    print("Saving deterministic AMR policies data to: ", filename)    
    df3.to_csv(filename, index=False)
    filename = output_dir + "/MO_eval.csv"
    print("Saving multi-objective eval loss to: ", filename)
    df4.to_csv(filename, index=False)

"""
    STEP 5: Plots
"""

if plot_figs or save_figs:

    import subprocess
    print("Calling plots.py")
    string_to_call = "python plots.py " + output_dir
    subprocess.call(string_to_call, shell=True)

    # print name of output_dir to file for plotting with slurm scritps
    file = open("output_dir4plots.txt","a")
    file.write("\n" + output_dir)
    file.close() 

print("\n\n\n ****** WARNING: hard coding initial parameters \n\n")
