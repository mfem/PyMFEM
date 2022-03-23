"""
    
    EXAMPLE 1-B: h-refinement policy for L-shaped domain problem with robust thresholds

"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from tensorflow.python.ops.gen_array_ops import lower_bound
from prob_envs.RobustPoisson import RobustThresholdPoisson
import numpy as np
import time
import seaborn as sns

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

train = True
nbatches = 250
minimum_budget_problem = False
save_figs = True

## Configuration for minimum budget problem
prob_config = {
    'mesh_name'             : 'l-shape-benchmark.mesh',
    'problem_type'          : 'lshaped',
    'num_unif_ref'          : 1,
    'order'                 : 2,
    'optimization_type'     : 'error_threshold', 
    'dof_threshold'         : 5e5,
    'error_threshold'       : 1e-4,
    'error_threshold_upper' : 5e-3,
    'error_threshold_lower' : 1.5e-4,
}

## Change to minimum error problem
if not minimum_budget_problem:
    prob_config['optimization_type'] = 'dof_threshold'
    prob_config['dof_threshold']     = 1e4
    prob_config['dof_threshold_lower'] = 5e3
    prob_config['dof_threshold_upper'] = 1.5e4

# Neural network configuration
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

"""
    STEP 2: Training
"""

## Train policy
ray.shutdown()
ray.init(ignore_reinit_error=True)
register_env("my_env", lambda config : RobustThresholdPoisson(**prob_config))
agent = ppo.PPOTrainer(env="my_env", config=config)
env = RobustThresholdPoisson(**prob_config)
if train:
    env.trainingmode = True
    for n in range(nbatches):
        print("training batch %d of %d batches" % (n+1,nbatches))
        result = agent.train()
        episode_score = -result["episode_reward_mean"]
        episode_length = result["episode_len_mean"]
        print ("Mean episode cost:   %.3f" % episode_score)
        print ("Mean episode length: %.3f" % episode_length)
    checkpoint_path = agent.save()
    print(checkpoint_path)
else:
    checkpoint_path = "/Users/keith10/ray_results/PPO_my_env_2022-01-07_21-22-29rwe7vrq4/checkpoint_000250/checkpoint-250"
    # checkpoint_path = "/Users/keith10/ray_results/PPO_my_env_2022-01-07_18-21-14a51g436k/checkpoint_000250/checkpoint-250"
    # checkpoint_path = "/Users/keith10/ray_results/PPO_my_env_2022-01-07_17-38-26wu1ex1wi/checkpoint_000250/checkpoint-250"

"""
    STEP 3: Validation
"""

error_threshold = prob_config['error_threshold']
dof_threshold = 1e4
# dof_threshold = prob_config['dof_threshold']
env.set_error_threshold(error_threshold)
env.set_dof_threshold(dof_threshold)

## Enact trained policy
agent.restore(checkpoint_path)
rlactions = []
obs = env.reset(random_threshold=False)
done = False
rlepisode_cost = 0
rlerrors = [env.global_error]
rldofs = [env.sum_of_dofs]
env.trainingmode = False
while not done:
    action = agent.compute_action(obs,explore=False)
    obs, reward, done, info = env.step(action)
    if not minimum_budget_problem and done:
       break
    rlactions.append(action)
    rlepisode_cost -= reward
    print("step = ", env.k)
    print("action = ", action.item())
    print("Num. Elems. = ", env.mesh.GetNE())
    print("episode cost = ", rlepisode_cost)
    rldofs.append(info['num_dofs'])
    rlerrors.append(info['global_error'])
    
    # time.sleep(0.05)
    # env.RenderMesh()

## Enact AMR policies
costs = []
actions = []
nth = 100
errors = []
dofs = []
for i in range(1, nth):
    print(i)
    action = np.array([i/nth])
    actions.append(action.item())
    obs = env.reset(random_threshold=False)
    done = False
    episode_cost_tmp = 0
    errors_tmp = [env.global_error]
    dofs_tmp = [env.sum_of_dofs]
    while not done:
        _, reward, done, info = env.step(action)
        if not minimum_budget_problem and done:
            break
        episode_cost_tmp -= reward
        errors_tmp.append(info['global_error'])
        dofs_tmp.append(info['num_dofs'])
    costs.append(episode_cost_tmp)
    errors.append(errors_tmp)
    dofs.append(dofs_tmp)
    

"""
    STEP 4: Plots
"""

## Plot training curve
root_path, _ = os.path.split(checkpoint_path)
root_path, _ = os.path.split(root_path)
csv_path = root_path + '/progress.csv'
df = pd.read_csv(csv_path)
cost = -df.episode_reward_mean.to_numpy()

plt.figure()
plt.plot(cost, color=palette_list[3], lw=2, label=r'Training curve')
ax1 = plt.gca()
ax1.set_xlabel("Epoch")
if minimum_budget_problem:
    ax1.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
else:
    ax1.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
ax1.legend()
plt.savefig('Example1b_fig1.pdf',format='pdf', bbox_inches='tight')

## Make letter-box plot
plt.figure()
ax2 = sns.boxenplot(y=costs, width=.6, color=palette_list[3], label='_nolegend_')
x2 = ax2.get_xlim()
plt.hlines(rlepisode_cost, x2[0], x2[1], lw=2, color='k', label=r'(AM)$^2$R policy cost')
y2 = ax2.get_ylim()
plt.fill_between(x2, np.floor(y2[0]), rlepisode_cost, color=palette_list[0], label=r'Apparent performance barrier')

ax2.set_xlabel(r'')
if minimum_budget_problem:
    ax2.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
else:
    ax2.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
lgd = ax2.legend()
letterbox_entry(lgd)
sns.despine(ax=ax2, bottom=True)
ax2.tick_params(bottom=False)
plt.tight_layout()
plt.savefig('Example1b_fig2.pdf',format='pdf', bbox_inches='tight')

## Plot theta vs. cost
plt.figure()
ax3 = plt.gca()
plt.plot(actions[9::10], costs[9::10], 'o', color=palette_list[3], label=r'AMR policies')
x3 = ax3.get_xlim()
plt.hlines(rlepisode_cost, x3[0], x3[1], lw=2, color='k', label=r'(AM)$^2$R policy')
y3 = ax3.get_ylim()
plt.fill_between(x3, np.floor(y3[0]), rlepisode_cost, color=palette_list[0], label=r'Apparent performance barrier')

ax3.set_xlabel(r'$\theta$ (constant) in AMR policy')
if minimum_budget_problem:
    ax3.set_ylabel(r'$\log_2(N_{\rm{dofs}})$')
else:
    ax3.set_ylabel(r'$\log_2(E_{\rm{Final}})$')
ax3.legend(loc='upper center')
plt.tight_layout()
plt.savefig('Example1b_fig3.pdf',format='pdf', bbox_inches='tight')

## Make convergence plots (1/2)
plt.figure()
ax4 = plt.gca()
alpha = 1.0
plt.loglog(dofs[9],errors[9],'-o',lw=1.3, color=palette_list[3], alpha=alpha, label=r'AMR policies')
# plt.loglog(dofs[9][-1],errors[9][-1], marker="o", markersize=10, color=palette_list[3], label='_nolegend_')
for k in range(19,len(errors),10):
    plt.loglog(dofs[k],errors[k],'-o',lw=1.3, color=palette_list[3], label='_nolegend_')
    # plt.loglog(dofs[k][-1],errors[k][-1], marker="o", markersize=10, color=palette_list[3], alpha=alpha, label='_nolegend_')
plt.loglog(rldofs,rlerrors,'-o',lw=1.3, color=palette_list[0], label=r'(AM)$^2$R policy')
# plt.loglog(rldofs[-1],rlerrors[-1], marker="o", markersize=10, color=palette_list[0], label='_nolegend_')
ax4.set_xlabel(r'Degrees of freedom')
ax4.set_ylabel(r'Relative error')
ax4.legend()
plt.savefig('Example1b_fig4.pdf',format='pdf', bbox_inches='tight')

## Make convergence plots (2/2)
cumdofs = []
cumrldofs = np.cumsum(rldofs)
for k in range(len(dofs)):
    cumdofs.append(np.cumsum(dofs[k]))
plt.figure()
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
plt.savefig('Example1b_fig5.pdf',format='pdf', bbox_inches='tight')

## Plot action vs. refinement step
plt.figure()
ax6 = plt.gca()
plt.plot(rlactions,'-o',lw=1.3, label=r'(AM)$^2$R policy')
ax6.set_xlabel(r'Refinement step')
ax6.set_ylabel(r'$\theta$ in trained (AM)$^2$R policy')
plt.savefig('Example1b_fig6.pdf',format='pdf', bbox_inches='tight')


plt.show()