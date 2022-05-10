"""
    
    EXAMPLE 2-C: hp-refinement policy for L-shaped domain problem with robust angles

"""

import os
from sqlite3 import DatabaseError
import matplotlib.pyplot as plt
import pandas as pd
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from tensorflow.python.ops.gen_array_ops import lower_bound
from prob_envs.RobustPoisson import hpPoisson
import numpy as np
import time
import seaborn as sns

from ray.tune.logger import UnifiedLogger
from datetime import datetime
import json 
import csv

def print_config(dir, prob_config = None, rl_config = None):
    if (prob_config is not None):
        with open(dir+"/prob_config.json", 'w') as f: 
            json.dump(prob_config,f)
            # for key, value in prob_config.items(): 
                # f.write('%s:%s\n' % (key, value))

    if (rl_config is not None):
        with open(dir+"/rl_config.json", 'w') as f: 
            json.dump(prob_config,f)

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
parser.add_argument('--marginals', dest='marginals_eval', default=False, action='store_true')
parser.add_argument('--angle', dest='angle_eval', type=float, default= np.pi / 2)
args = parser.parse_args()
print("\n Parsed options = ", args, "\n")
train=args.train
eval=args.eval
save_data=args.savedata
plot_figs=args.plotfigs
save_figs=args.savefigs
angle_abbrv = "{:.2f}".format(np.round(args.angle_eval,2)) # to fix filename length

restore_policy = False
nbatches = 100 #250
minimum_budget_problem = False  # mininum budget == error_threshold == minimize dofs

## Configuration for minimum budget problem
prob_config = {
    'mesh_name'             :  'circle_3_4.mesh', # 'l-shape-benchmark.mesh', 'star.mesh', 'circle_3_4.mesh'
    'problem_type'          :  'lshaped', #'noneoftheabove', 
    'num_unif_ref'          : 1,
    'order'                 : 2,
    'optimization_type'     : 'error_threshold', 
    'dof_threshold'         : 5e5,
    'error_threshold'       : 1e-4,
    'angle_lower'           : np.pi * 0.1, # np.pi * 0.25,
    'angle_upper'           : np.pi * 1.9, # np.pi * 0.75,
    'num_batches'           : nbatches
}

## Change to minimum error problem
if not minimum_budget_problem:
    prob_config['optimization_type'] = 'dof_threshold'
    prob_config['dof_threshold']     = 1e4



## Neural network configuration
model_config = {
    "fcnet_hiddens"    : [128, 128],
    "fcnet_activation" : "swish",
}

## rllib parameters
config = ppo.DEFAULT_CONFIG.copy()
# config['batch_mode'] = 'truncate_episodes'
config['batch_mode'] = 'complete_episodes'
config['sgd_minibatch_size'] = 100
config['rollout_fragment_length'] = 50
config['num_workers'] = 10
config['train_batch_size'] = config['rollout_fragment_length'] * config['num_workers']
config['num_gpus'] = 0
config['gamma'] = 1.0
config['lr'] = 1e-4
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
    # temp_path = 'Example2c_2022-05-06_19-33-30/' # dof thresh 1e4; L-shape mesh; pi/4 - 3pi/4 angle training; num_batch 250
    # temp_path = 'Example2c_2022-05-09_14-23-16/' # dof thresh 1e4; pacman mesh; 0.5 - 1.5 pi angle training; num_batch 100
    temp_path = 'Example2c_2022-05-09_17-32-28/' # dof thresh 1e4; pacman mesh; 0.1 - 1.9 pi angle training; num_batch 100
    

    
    checkpoint_dir = log_dir + temp_path
    chkpt_file=checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num)
    output_dir = output_dir_ + temp_path
else:
    timestr = datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    temp_path = 'Example2c_' + timestr
    checkpoint_dir = log_dir + temp_path
    output_dir = output_dir_ + temp_path

## Train policy
ray.shutdown()
ray.init(ignore_reinit_error=True)
register_env("my_env", lambda config : hpPoisson(**prob_config))
trainer = ppo.PPOTrainer(env="my_env", config=config,
                    logger_creator=custom_log_creator(checkpoint_dir))
env = hpPoisson(**prob_config)

# obs = env.reset()
# env.render()
# action = np.array([0.1, 0.1])
# obs, reward, done, info = env.step(action)
# env.render()
# exit()

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
    checkpoint_path = trainer.save()
    print(checkpoint_path)
if eval and not train:
    # temp_path = 'Example2c_2022-05-06_19-33-30/' # dof thresh 1e4; L-shape mesh; pi/4 - 3pi/4 angle training; num_batch 250
    # temp_path = 'Example2c_2022-05-09_14-23-16/' # dof thresh 1e4; pacman mesh; 0.5 - 1.5 pi angle training; num_batch 100
    temp_path = 'Example2c_2022-05-09_17-32-28/' # dof thresh 1e4; pacman mesh; 0.1 - 1.9 pi angle training; num_batch 100
    chkpt_num = nbatches
    checkpoint_dir = log_dir + temp_path
    # checkpoint_path=checkpoint_dir+'/checkpoint_0000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt < 100
    checkpoint_path=checkpoint_dir+'/checkpoint_000'+str(chkpt_num)+'/checkpoint-'+str(chkpt_num) # if checkpt > 99 and <1000
    output_dir = output_dir_ + temp_path

if train:
    print_config(checkpoint_dir,prob_config=prob_config, rl_config=rl_config)

"""
    STEP 3: Validation
"""

if eval and not args.marginals_eval:

    if prob_config['mesh_name'] == 'l-shape-benchmark.mesh' or prob_config['mesh_name'] == 'circle_3_4.mesh': 
        env.set_angle(args.angle_eval) # note: set_angle(...) redefines self.initial_mesh
        print("*** Set angle for eval to  ", args.angle_eval)

    ## Enact trained policy
    trainer.restore(checkpoint_path)
    cols = ['theta', 'rho']
    rlactions = pd.DataFrame(columns=cols)
    obs = env.reset(random_angle=False)
    done = False
    rlepisode_cost = 0
    rlerrors = [env.global_error]
    rldofs = [env.sum_of_dofs]
    env.trainingmode = False
    num_steps_of_RL_policy = 0


    while not done:
        action = trainer.compute_single_action(obs,explore=False)
        obs, reward, done, info = env.step(action)
        if not minimum_budget_problem and done:
            break
        rlactions = rlactions.append({'theta': action[0], 'rho':action[1]}, ignore_index=True)
        rlepisode_cost -= reward
        print("step = ", env.k)
        print("theta action = ", action[0].item())
        print("rho   action = ", action[1].item())
        print("Num. Elems. = ", env.mesh.GetNE())
        print("episode cost = ", rlepisode_cost)
        rldofs.append(info['num_dofs'])
        rlerrors.append(info['global_error'])
        # print()
        # np.set_printoptions(precision=16)
        # print(env.k)
        # print(np.sort(env.errors))
        # env.RenderHPmesh()
        # if env.k > 6:
        #     exit()

    # print("*** done with RL eval - exiting")
    # exit()
    
    # env.RenderMesh()
    # env.RenderHPmesh()
    # print("\nRendering RL policies and exiting\n")
    # exit()
    
    ## Enact AMR policies, using "expert" strategy
    costs = []
    actions = []
    nth = 100
    errors = []
    dofs = []
    for i in range(1, nth):
        print(i)
        action = np.array([i/nth]) # note 1D action space
        actions.append(action.item())
        print("action = ", action.item())
        obs = env.reset(random_angle=False)
        done = False
        episode_cost_tmp = 0
        errors_tmp = [env.global_error]
        dofs_tmp = [env.sum_of_dofs]
        max_steps   = 40 # or:  num_steps_of_RL_policy
        steps_taken = 0
        while not done:
            _, reward, done, info = env.expertStep(action)
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

    # ## Enact AMR policies, using "fixed two parameter sweep" strategy
    tp_costs = []
    tp_nth = 10
    # tp_actions = np.zeros(((tp_nth-1)**2,2)) # exclude 0.0, 1.0 as actions
    # tp_actions = np.zeros(((tp_nth+1)**2,2)) # include 0.0, 1.0 as actions
    tp_actions = np.zeros(((tp_nth)**2,2)) # include 0.0 but not 1.0 as actions
    tp_errors = []
    tp_dofs = []
    index_count = 0

    # for theta in range(1, tp_nth):      # exclude 0.0, 1.0 as actions
    #     for rho in range(1, tp_nth):    # exclude 0.0, 1.0 as actions
    # for theta in range(0, tp_nth+1):    # include 0.0, 1.0 as actions
    #     for rho in range(0, tp_nth+1):  # include 0.0, 1.0 as actions
    for theta in range(0, tp_nth):        # include 0.0 but not 1.0 as actions
        for rho in range(0, tp_nth):      # include 0.0 but not 1.0 as actions
            tp_actions[index_count] = np.array([theta/tp_nth, rho/tp_nth]) # note 1D action space
            if theta/tp_nth == 1 and rho/tp_nth == 1: # avoid some linear algerbra error if action is [1,1]
                tp_actions[index_count] = np.array([0.99, 0.99])
            action = tp_actions[index_count]
            print("action = ", tp_actions[index_count])
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
            print('episode cost = ', episode_cost_tmp)
            if theta == 5 and rho == 5:
                env.RenderMesh()
                env.RenderHPmesh()
                print("\nRendering two parmaeter policy ", tp_actions[index_count-1], "\n")

if args.marginals_eval:
    print("Creating data for marginals")
    mkdir_p(output_dir+"marginals/")

    # angle_vals = np.pi* np.linspace(1/4,3/4,3)    # 3 angles, for debugging
    # angle_vals = np.pi* np.linspace(1/4,3/4,100)  # 100 angle test
    # angle_vals = np.pi* np.linspace(1/2,1/2,1)    # only do pi/2
    # angle_vals = np.pi* np.linspace(0,0,1)        # only do 0
    # angle_vals = np.pi* np.linspace(3/4,95/100,10)    
    # # angle_vals = np.pi* np.linspace(5/100,1/4,10)    
    # angle_vals = np.pi* np.linspace(5/100,95/100,20)  
    # angle_vals = np.pi* np.linspace(0,1,21) # from 0 to pi in increments of 0.05 * pi
    # angle_vals = np.pi* np.linspace(1/4,3/4,21)  # pi/4 to 3pi/4 only
    angle_vals = np.pi* np.linspace(0.1, 1.9, 21)  # .1 to 1.9 pi
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
    filename = output_dir+"marginals/marginals_rl.csv"
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


    # import time
    # job_time_id = int(time.time()) 
    # print("Job time id = ", job_time_id)

    #### tpp marginals data df and saving
    df_tpp = pd.DataFrame(rows, columns=headers)
    filename = output_dir+"marginals/marginals_tpp.csv"
    print("Saving marginals data to: ", filename)    
    df_tpp.to_csv(filename, index=False)
    print("Exiting")
    exit()




  
"""
    STEP 4: Save Data
"""

if (save_data or save_figs):
    mkdir_p(output_dir)
    print()
    print_config(output_dir,prob_config=prob_config, rl_config=rl_config)
    print()


if save_data:
    import json
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
    filename = output_dir+"/rl_data_angle_" + angle_abbrv + ".csv"
    print("\nSaving RL deployed policy data to: ", filename)  
    df2.to_csv(filename, index=False)

    #### expert policy df
    df3 = pd.DataFrame({'actions':actions,'costs':costs,'errors':errors,'dofs':dofs})
    filename = output_dir+"/deterministic_amr_data_angle_" + angle_abbrv + ".csv"
    print("\nSaving deterministic AMR policies data to: ", filename)    
    df3.to_csv(filename, index=False)

    #### two param policy df
    df4 = pd.DataFrame({'theta':tp_actions[:,0], 'rho':tp_actions[:,1],'costs':tp_costs,'errors':tp_errors,'dofs':tp_dofs})
    filename = output_dir+"/two_param_amr_data_angle_" + angle_abbrv + ".csv"
    print("\nSaving two parameter AMR policies data to: ", filename)    
    df4.to_csv(filename, index=False)

"""
    STEP 5: Plots
"""

if plot_figs or save_figs:

    import subprocess
    print("Calling plots.py")
    string_to_call = "python plots.py " + output_dir + " --angle " + angle_abbrv
    subprocess.call(string_to_call, shell=True)