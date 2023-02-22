"""
    
    EXAMPLE 2-C multi-objective: hp-refinement policy for L-shaped domain problem with robust angles with multi-objective cost function

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from prob_envs.MultiObjectivePoisson import hp_Angle_MultiObjPoisson
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

parser.add_argument('--marginals', dest='marginals_eval', default=False, action='store_true')
parser.add_argument('--angle', dest='angle_eval', type=float, default= np.pi / 2)
parser.add_argument('--meshnum', dest='command_line_mesh', type=int, default=0)

parser.add_argument('--observe_budget', default = True, action='store_true')
parser.add_argument('--no_observe_budget', dest='observe_budget', action='store_false')


args = parser.parse_args()
print("Parsed options = ", args)

train       = args.train
eval        = args.eval
save_data   = args.savedata
plot_figs   = args.plotfigs
save_figs   = args.savefigs
angle_abbrv = "{:.2f}".format(np.round(args.angle_eval,2)) # to keep filename length short

restore_policy = False
nbatches       = 50

## Configuration for multi objective problem
prob_config = {
    'mesh_name'         : 'l-shape-benchmark.mesh', #'star.mesh',
    'problem_type'      : 'lshaped', #'default',
    'num_unif_ref'      : 1,
    'order'             : 2,
    'optimization_type' : 'multi_objective', 
    'cost_function'     : 'ADF', # ADF for annealing desirability function
    'dof_threshold'     : 1e6,
    'angle_lower'       : np.pi * 0.1, # np.pi * 0.25,
    'angle_upper'       : np.pi * 0.9, # np.pi * 0.75,

    'observe_budget'    : args.observe_budget,

    'num_iterations'    : 50,
    'num_batches'       : nbatches,
    'tau_min'           : np.log2(5e-4), #np.log2(1e-2), # np.log2(1e-3), #np.log2(1e-4),
    'tau_max'           : np.log2(5e-2),  #np.log2(5e-2),
    'M_warm'            : 5, #50 number of batches in warming phase
    'M_anneal'          : 5,  #25, # number of batches per tau in annealing phase
    'N_anneal'          : 1,  #10,  # number of target errors to train on (not counting intitial target)
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
config['batch_mode']            = 'truncate_episodes'
# config['batch_mode']              = 'complete_episodes'
config['sgd_minibatch_size']      = 100
config['rollout_fragment_length'] = 50
config['num_workers']             = 10
config['train_batch_size']        = config['rollout_fragment_length'] * config['num_workers']
config['num_gpus']                = 0
config['gamma']                   = 1.0
config['lr']                      = 1e-4 
config['seed']                    = 4000
config['model']                   = model_config


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

# Define other meshes
# if command_line_mesh == 0, prob_config is not altered by this block
if args.command_line_mesh == 1:
    prob_config['mesh_name'] = 'circle_3_4.mesh'
    prob_config['num_unif_ref'] = 0
if args.command_line_mesh == 2:
    prob_config['mesh_name'] = 'star.mesh'
    prob_config['problem_type'] = 'default'
    prob_config['num_unif_ref'] = 0
if args.command_line_mesh == 3:
    prob_config['mesh_name'] = 'staircase.mesh'
    prob_config['problem_type'] = 'default'
    prob_config['num_unif_ref'] = 0
if args.command_line_mesh == 4:
    prob_config['mesh_name'] = 'staircase_tri2.mesh'
    prob_config['problem_type'] = 'default'
    prob_config['num_unif_ref'] = 0
if args.command_line_mesh == 5:
    prob_config['mesh_name'] = 'staircase_tri.mesh'
    prob_config['problem_type'] = 'default'
    prob_config['num_unif_ref'] = 0
if args.command_line_mesh == 6:
    prob_config['mesh_name'] = 'fichera.mesh'
    prob_config['problem_type']  = 'default'
    prob_config['num_unif_ref']  = 0
    prob_config['dof_threshold'] = 5e6  # note higher dof threshold for Fichera mesh
    
if prob_config['mesh_name'] == 'l-shape-benchmark.mesh':
    mesh_abbrv  = 'lshape'
elif prob_config['mesh_name'] == 'circle_3_4.mesh':
    mesh_abbrv  = 'pacman'
elif prob_config['mesh_name'] == 'star.mesh':
    mesh_abbrv  = 'star'
    angle_abbrv = 'nan'
elif prob_config['mesh_name'] == 'fichera.mesh':
    mesh_abbrv  = 'fichera'
    angle_abbrv = 'nan'
elif prob_config['mesh_name'] == 'staircase.mesh':
    mesh_abbrv = 'stair'
    angle_abbrv = 'nan'
elif prob_config['mesh_name'] == 'staircase_tri.mesh':
    mesh_abbrv = 'stairtri1'
    angle_abbrv = 'nan'
elif prob_config['mesh_name'] == 'staircase_tri2.mesh':
    mesh_abbrv = 'stairtri2'
    angle_abbrv = 'nan'
elif prob_config['mesh_name'] == 'square-disc.mesh':
    mesh_abbrv = 'squaredisc'
    angle_abbrv = 'nan'
else:
    mesh_abbrv = prob_config['mesh_name']
    angle_abbrv = 'nan'


"""
    STEP 2: Training
"""

homepath = os.path.expanduser("~")
log_dir = os.getcwd() + '/logs/'
output_dir_ = os.getcwd() + '/output/'

if (restore_policy):
    chkpt_num = nbatches
    # set the path of the checkpoint
    temp_path = ''
    checkpoint_dir = log_dir + temp_path
    chkpt_file=checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num)
    output_dir = output_dir_ + temp_path
else:
    timestr = datetime.today().strftime("%Y-%m-%d_%H-%M-%S")

    if prob_config['optimization_type'] == 'multi_objective':
        if ADF_bool:
            temp_path = 'Example2c_ADF_' + timestr
        
    else:
        temp_path = 'Example2c_' + timestr

    checkpoint_dir = log_dir + temp_path
    output_dir = output_dir_ + temp_path   

## Train policy
ray.shutdown()
ray.init(ignore_reinit_error=True)


# initialize parameter class
ADF_params = ADF_Parameters.options(name = "parameters").remote()

register_env("my_env", lambda config : hp_Angle_MultiObjPoisson(**prob_config))
trainer = ppo.PPOTrainer(env="my_env", config=config, 
                       logger_creator=custom_log_creator(checkpoint_dir))
env = hp_Angle_MultiObjPoisson(**prob_config)
print("env = ", env)
print("prob config = ", prob_config)



if (restore_policy):
   trainer.restore(chkpt_file)

# initialize tau and delta in ADF_Parameter class
if ADF_bool:
    tau_min = prob_config['tau_min']
    tau_max = prob_config['tau_max']
    
    if N_anneal > 0:
        tau_step = (tau_max - tau_min)/N_anneal;
    else:
        tau_step = 0; 

    delta_warm = 1
    delta_anneal = 1

    print("assiging params: ", tau_max, tau_step, delta_warm, delta_anneal)
    ADF_params.set_initial_parameters.remote(tau_max, tau_step, delta_warm, delta_anneal)

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
        print("config: ")
        print(config)
        print("rl_config:")
        print(rl_config)
        result = trainer.train()
        episode_score = -result["episode_reward_mean"]
        episode_length = result["episode_len_mean"]
        print ("Mean episode cost:   %.3f" % episode_score)

        checkpoint_path = trainer.save(checkpoint_dir)
        print(checkpoint_path)

if eval and not train:
    MO_eval_loss = []
    if prob_config['optimization_type'] == 'multi_objective':
        temp_path = ''
    else:
        print("*** need to set path for eval *** exiting")
        exit()


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

# determine which errors/taus to evaluate on (tau is the 2log of the target error)

if eval:
    env.set_angle(args.angle_eval) # note: set_angle(...) redefines self.initial_mesh
    print("*** Set angle for eval to  ", args.angle_eval)

    if N_anneal > 0:
        params_to_eval = np.array([tau_max - i * tau_step for i in range(N_anneal + 1)])
    else:
        params_to_eval = np.array([tau_max])
    #params_to_eval = np.array([np.log2(0.01), np.log2(0.005), np.log2(0.001)])
    params_to_eval = np.array([params_to_eval[5], params_to_eval[8], params_to_eval[10]])
    print("\n *** params to eval =", params_to_eval, "\n")

else:
    params_to_eval = []

if eval:
    # create file location for saving evaluation results
    mkdir_p(output_dir)
    file_name = "/ADF_policy_data_2c.txt"
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

    adf_vs_fixed_df = pd.DataFrame(columns=['target', 'theta', 'rho', 'steps', 'dofs', 'error', 'why'])

    save_every_step = True;
    for param in params_to_eval:
        ## Enact trained policy
        trainer.restore(checkpoint_path) 
        rlactions = []

        # set tau
        ADF_params.set_tau.remote(param)
        tau_str = str(param).replace('.','_')
        obs = env.reset(random_angle = False)
        
        done = False
        rlepisode_cost = 0
        rlerrors = [env.global_error]
        rldofs = [env.sum_of_dofs]
        env.trainingmode = False
        quit_reason = 'notdone'
        print("\n\n")

        while not done:
            # action = trainer.compute_single_action(obs, explore = False)
            action = trainer.compute_single_action(obs, explore=True)
            obs, reward, done, info = env.step(action)
            if prob_config['optimization_type'] == 'dof_threshold' and done:
                break
            rlactions.append({'theta' : action[0], 'rho' : action[1]}, ignore_index = True) 
            rlepisode_cost -= reward
            print("RL step =", env.k, 
                  "theta =", np.round(float(action[0].item()),3),
                  "rho =", np.round(float(action[1].item()),3),
                  "dof=", info['num_dofs'], 
                  # dofs at this step = env.mesh.GetNE(),
                  "err=", np.round(info['global_error'],6),
                  "cum cost=", rlepisode_cost)
            rldofs.append(info['num_dofs']) 
            rlerrors.append(info['global_error'])
            quit_reason = info['why']

            if save_every_step == True and args.savemesh:
                print("step {}".format(env.k))
                env.render()
                env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + prob_config['problem_type'] + "_tau_" + tau_str + '_step_' + str(env.k) + '.mesh')

        save_every_step = False; # set this to false after saving one set of refinement meshes

        print("")


        if args.savemesh:

            env.render()
            env.mesh.Save(output_dir+"/meshes_and_gfs/" + 'rl_mesh_' + prob_config['problem_type'] + "_tau_" + tau_str + '.mesh')

        # save final errors in file for each tau
        cum_rldofs = np.cumsum(rldofs)

        file_string = str(ray.get(ADF_params.get_tau.remote())) + ", " + str(cum_rldofs[-1]) + ", " + str(rlerrors[-1]) + ", " + str(env.k) + "\n"
        file = open(file_location, "a")
        file.write(file_string)
        file.close()

        # save RL data to dataframe
        adf_vs_fixed_df = adf_vs_fixed_df.append(
            {'target': np.round(np.power(2,param),8), 'theta': -1, 'rho' : -1, 'steps': env.k, 'dofs': env.sum_of_dofs, 'error' : rlerrors[-1], 'why': quit_reason},
            ignore_index=True)
        
        dftgt = pd.DataFrame({'rlactions':rlactions})
        current_target_error = np.round(np.power(2,param),8)
        filename = output_dir+"/rl_data_tgt_" + str(current_target_error) + ".csv"
        print("Saving RL actions for error target ", current_target_error, "to: ", filename)    
        dftgt.to_csv(filename, index=False)
        
                
        #############
        # Evaluate "expert AMR policies, if L-shaped problem type
        ########### 
        print("\n\n *** skipping expert AMR policies - hard coded temporarily; remove `and false` in two places below *** \n\n")

        if prob_config['problem_type'] == 'lshaped' and False:
            ## Enact AMR policies
            costs = []
            actions = []
            nth = 100 # number of fixed policies to try
            errors = []
            dofs = []
            for i in range(1, nth):
                quit_reason = 'notdone'
                action = np.array([i/nth])
                actions.append(action.item())
                # print("action = ", action.item())
             
                obs = env.reset(random_angle = False)
                done = False
                episode_cost_tmp = 0
                errors_tmp = [env.global_error]
                dofs_tmp = [env.sum_of_dofs]
                
                while not done:
                    _, reward, done, info = env.expertStep(action)
                    episode_cost_tmp -= reward
                    errors_tmp.append(info['global_error'])
                    dofs_tmp.append(info['num_dofs'])
                    quit_reason = info['why']
                costs.append(episode_cost_tmp)
                errors.append(errors_tmp)
                dofs.append(dofs_tmp)
                # print('episode cost = ', episode_cost_tmp)
                adf_vs_fixed_df = adf_vs_fixed_df.append(
                    {'target': np.round(np.power(2,param),8), 
                    'theta'  : action[0],
                    'rho'    : action[1],
                    'steps'  : env.k,
                    'dofs'   : env.sum_of_dofs,
                    'error'  : errors[-1],
                    'why'    : quit_reason},
                    ignore_index=True)
                print("expert action ", action.item(), " had episode cost = :", episode_cost_tmp)


        fixed_parameter_eval = False;
        print("\n *** skipping fixed two parameter policies - hard coded temporarily.*** \n")
        if fixed_parameter_eval:
            #############
            # Enact fixed "two parameter" AMR policies
            ########### 
            tp_costs = []
            tp_nth = 10
            # tp_actions = np.zeros(((tp_nth-1)**2,2)) # exclude 0.0, 1.0 as actions
            # tp_actions = np.zeros(((tp_nth+1)**2,2)) # include 0.0, 1.0 as actions
            # tp_actions = np.zeros(((tp_nth)**2,2)) # include 0.0 but not 1.0 as actions
            tp_actions = np.zeros((1,2)) # include only theta = [one fixed value], rho = [one fixed value]
            tp_errors = []
            tp_dofs = []
            index_count = 0

            # for theta in range(1, tp_nth):      # exclude 0.0, 1.0 as actions
            #     for rho in range(1, tp_nth):    # exclude 0.0, 1.0 as actions
            # for theta in range(0, tp_nth+1):    # include 0.0, 1.0 as actions
            #     for rho in range(0, tp_nth+1):  # include 0.0, 1.0 as actions
            # for theta in range(0, tp_nth):        # include 0.0 but not 1.0 as actions
            #     for rho in range(0, tp_nth):      # include 0.0 but not 1.0 as actions
            # for theta in range(6, 7):        # include only theta = 0.6, rho = 0.5
            #     for rho in range(5, 6):      # include only theta = 0.6, rho = 0.5
            # for theta in range(6, 7):        # include only theta = 0.6, rho = 0.4
            #     for rho in range(4, 5):      # include only theta = 0.6, rho = 0.4
            for theta in range(6, 7):        # include only theta = 0.6, rho = 0.3
                for rho in range(3, 4):      # include only theta = 0.6, rho = 0.3
                    tp_actions[index_count] = np.array([theta/tp_nth, rho/tp_nth]) # note 1D action space
                    if theta/tp_nth == 1 and rho/tp_nth == 1: # avoid some linear algerbra error if action is [1,1]
                        tp_actions[index_count] = np.array([0.99, 0.99])
                    action = tp_actions[index_count]
                    print("two param action = ", tp_actions[index_count])
                    index_count += 1
                    obs = env.reset(random_angle=False)
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
          

    if args.marginals_eval:
        print("Creating data for marginals")
        mkdir_p(output_dir+"/marginals/")

        # angle_vals = np.pi* np.linspace(1/4,3/4,3)    # 3 angles, for debugging
        # angle_vals = np.pi* np.linspace(1/4,3/4,100)  # 100 angle test
        # angle_vals = np.pi* np.linspace(1/2,1/2,1)    # only do pi/2
        # angle_vals = np.pi* np.linspace(0,0,1)        # only do 0
        # angle_vals = np.pi* np.linspace(3/4,95/100,10)    
        # # angle_vals = np.pi* np.linspace(5/100,1/4,10)    
        # angle_vals = np.pi* np.linspace(5/100,95/100,20)  
        # angle_vals = np.pi* np.linspace(0,1,21) # from 0 to pi in increments of 0.05 * pi
        # angle_vals = np.pi* np.linspace(1/4,3/4,21)  # pi/4 to 3pi/4 only
        # angle_vals = np.pi* np.linspace(0.1, 1.9, 21)  # .1 to 1.9 pi
        angle_vals = np.pi* np.linspace(0.1, 0.9, 21)  # .1 to .9 pi
        print("Evaluating marginals with angle values = ", angle_vals)

       
        ##########################################
        # single RL policy marginals data creation 
        ##########################################
        
        headers = ['theta', 'rho', 'angle', 'step', 'num elts', 'num dofs', 'sum dofs', 'error est', 'cost']#, 'L2 Error', 'H1 Error']
        rows = []
     
        for angle in angle_vals:
            trainer.restore(checkpoint_path)
            env.set_angle(angle)
            print("*** Set angle for RL eval to  ", angle)
            obs = env.reset(random_angle=False)
            env.trainingmode = False
            done = False

            rlepisode_cost = 0.0
            steps_taken = 0
            max_steps = 40

            while not done:
                action = trainer.compute_single_action(obs,explore=False)
                rows.append([action[0], action[1], angle, env.k, env.mesh.GetNE(), env.fespace.GetTrueVSize(), 
                                            env.sum_of_dofs, env.global_error, rlepisode_cost])
                obs, reward, done, info = env.step(action)
                if not minimum_budget_problem and done:
                    break
                if steps_taken > max_steps:
                    print("*** BREAKING EARLY - fixed action exceeded max step threshold of ", max_steps, "steps.")
                    break
                else:
                    steps_taken += 1
                rlepisode_cost -= reward
                # print("step     = ", env.k)
                # print("action   = [", action[0], ", ", action[1], "]")
                # print("num elts = ", env.mesh.GetNE())
                # print("ep. cost = ", rlepisode_cost)
            
        
        #### rl marginals data df and saving
        df_rl = pd.DataFrame(rows, columns=headers)
        filename = output_dir+"/marginals/marginals_rl.csv"
        print("Saving RL marginals data to: ", filename, "\n")    
        df_rl.to_csv(filename, index=False)

        ##########################################
        # two param policy marginals data creation 
        ##########################################
       
        samples_per_param = 10

        headers = ['theta', 'rho', 'angle', 'step', 'num elts', 'num dofs', 'sum dofs', 'error est', 'cost']#, 'L2 Error', 'H1 Error']
        # headers = ['theta', 'rho',         'N', 'DoFs', 'Total_DoFs', 'Error_Estimate', 'Cost']#, 'L2_Error', 'H1_Error']
        rows = []
        for j in range(samples_per_param):
            print("... working on actions [", j/samples_per_param, ", :] ...")
            for k in range(samples_per_param):
                for angle in angle_vals:
                    env.set_angle(angle)
                    # print("*** Set angle for eval to  ", angle)
                    obs = env.reset(random_angle=False)
                    done = False
                    theta = j / samples_per_param
                    rho = k / samples_per_param
                    action = np.array([theta, rho])
                    # print("Collecting data for action ", action, " with angle ", angle)
                    episode_cost = 0.0
                    rows.append([action[0].item(), action[1].item(), angle, env.k, env.mesh.GetNE(), env.fespace.GetTrueVSize(), 
                                            env.sum_of_dofs, env.global_error, episode_cost])
                    steps_taken = 0
                    max_steps = 40
                    while not done:
                        _, reward, done, info = env.step(action)
                        episode_cost -= reward 
                        rows.append([action[0].item(), action[1].item(), angle, env.k, env.mesh.GetNE(), env.fespace.GetTrueVSize(), 
                                            env.sum_of_dofs, env.global_error, episode_cost])
                        if not minimum_budget_problem and done:
                            break
                        if steps_taken > max_steps:
                            print("*** BREAKING EARLY - fixed action exceeded max step threshold of ", max_steps, "steps.")
                            break
                        else:
                            steps_taken += 1
                            

                    # if options.save_mesh:
                    #     gfname = 'gf_theta_' + str(theta) + '_rho_' + str(rho) + '_angle_' + str(np.round(angle,5)) + '.gf'
                    #     env.RenderHPmesh(gfname=gfname)
                    #     env.mesh.Save('mesh_step_theta_' + str(theta) + '_rho_' + str(rho) + '_angle_' + str(np.round(angle,5)) + '.mesh')
                    #     # sol_sock = mfem.socketstream("localhost", 19916)
                    #     # sol_sock.precision(8)
                    #     # sol_sock.send_solution(self.mesh, orders)
                    #     # title = "step " + str(self.k)
                    #     # sol_sock.send_text('keys ARjlmp*******' + " window_title '" + title)
                    #     # sol_sock.send_text("valuerange 1.0 8.0 \n")
                    #     # sol_sock.send_text('keys S')

        #### tpp marginals data df and saving
        df_tpp = pd.DataFrame(rows, columns=headers)
        filename = output_dir+"/marginals/marginals_tpp.csv"
        print("Saving two parameter poilcy marginals data to: ", filename)    
        df_tpp.to_csv(filename, index=False)

        ##########################################
        # expert policy marginals data creation 
        ##########################################
       
        
        if prob_config['problem_type'] == 'lshaped':
            samples_per_param = 10

            headers = ['theta', 'angle', 'step', 'num elts', 'num dofs', 'sum dofs', 'error est', 'cost']#, 'L2 Error', 'H1 Error']
            # headers = ['theta', 'rho',         'N', 'DoFs', 'Total_DoFs', 'Error_Estimate', 'Cost']#, 'L2_Error', 'H1_Error']
            rows = []
            for j in range(samples_per_param):
                print("... working on expert action ", j/samples_per_param, " ...")
                for angle in angle_vals:
                    env.set_angle(angle)
                    obs = env.reset(random_angle=False)
                    done = False
                    theta = j / samples_per_param
                    action = np.array([theta]) # note 1D action space
                    episode_cost = 0.0
                    rows.append([action[0].item(), angle, env.k, env.mesh.GetNE(), env.fespace.GetTrueVSize(), 
                                            env.sum_of_dofs, env.global_error, episode_cost])
                    steps_taken = 0
                    max_steps = 40
                    while not done:
                        _, reward, done, info = env.expertStep(action)
                        episode_cost -= reward 
                        rows.append([action[0].item(), angle, env.k, env.mesh.GetNE(), env.fespace.GetTrueVSize(), 
                                                env.sum_of_dofs, env.global_error, episode_cost])
                        if not minimum_budget_problem and done:
                            break
                        if steps_taken > max_steps:
                            print("*** BREAKING EARLY - fixed action exceeded max step threshold of ", max_steps, "steps.")
                            break
                        else:
                            steps_taken += 1


            #### expert marginals data df and saving
            df_exp = pd.DataFrame(rows, columns=headers)
            filename = output_dir+"/marginals/marginals_exp.csv"
            print("Saving expert poilcy marginals data to: ", filename)    
            df_exp.to_csv(filename, index=False)

            print("Exiting")
            exit()


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
    df2 = pd.DataFrame({'theta':rlactions.iloc[:,0], 'rho' : rlactions.iloc[:,1],'rldofs':rldofs,'rlerrors':rlerrors})
    df3 = pd.DataFrame({'actions':actions,'costs':costs,'errors':errors,'dofs':dofs})
    df4 = pd.DataFrame({'MO_eval_loss':MO_eval_loss})

    filename = output_dir+"/training_data.csv"
    print("Saving training data to: ", filename)    
    df1.to_csv(filename, index=False)

    filename = output_dir+"/rl_data_" + mesh_abbrv + "_angle_" + angle_abbrv + ".csv"
    print("Saving RL deployed policy data to: ", filename)    
    df2.to_csv(filename, index=False)

    filename = output_dir+"/expert_data_" + mesh_abbrv + "_angle_" + angle_abbrv + ".csv"
    print("Saving deterministic AMR policies data to: ", filename)    
    df3.to_csv(filename, index=False)

    filename = output_dir + "/tpp_data_" + mesh_abbrv + "_angle_" + angle_abbrv + ".csv"
    print("Saving multi-objective eval loss to: ", filename)
    df4.to_csv(filename, index=False)

    filename = output_dir + "/adf_vs_fxd.csv"
    print("Saving ADF vs fixed data to:", filename)
    adf_vs_fixed_df.to_csv(filename, index=False)

"""
    STEP 5: Plots
"""

if plot_figs or save_figs:

    import subprocess
    print("Calling plots.py")
    string_to_call = "python plots.py " + output_dir
    subprocess.call(string_to_call, shell=True)

    '''
    # print name of output_dir to file for plotting with slurm scritps
    file = open("output_dir4plots.txt","a")
    file.write("\n" + output_dir)
    file.close() 
    '''
