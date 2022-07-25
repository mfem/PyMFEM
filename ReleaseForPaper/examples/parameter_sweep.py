import os
import matplotlib.pyplot as plt
import ray
import ray.rllib.agents.ppo as ppo
from ray.tune.registry import register_env
from prob_envs.MultiObjectivePoisson import MultiObjPoisson
import numpy as np
import time
import seaborn as sns

from ray.tune.logger import UnifiedLogger
from datetime import datetime
import json 

## Parameters
output_file_name 	= "fixed_theta_1a_AMR_data.txt"
num_iterations 		= 10

# values of alpha and theta to test
alphas = 0.1*np.array(range(0, 11, 1))
thetas = 0.1*np.array(range(1, 10, 1))

# define environment
prob_config = {
    'mesh_name'         : 'l-shape-benchmark.mesh',
    'problem_type'      : 'lshaped',
    'num_unif_ref'      : 1,
    'order'             : 2,
    'optimization_type' : 'multi_objective', 
    'dof_threshold'     : 1e5,
    'alpha'             : 0.5,
    'observe_alpha'     : False,
    'observe_budget'    : True,
    'num_iterations'    : 10
}
ray.shutdown()
ray.init(ignore_reinit_error=True)
register_env("my_env", lambda config : MultiObjPoisson(**prob_config))
#trainer = ppo.PPOTrainer(env="my_env", config=config, logger_creator=custom_log_creator(checkpoint_dir))
env = MultiObjPoisson(**prob_config)
env.trainingmode = False

# open output file
output_file = open(output_file_name, "a")
output_file.write("alpha, theta, cum_dofs, error\n")

# initialize variables
costs = []
actions = []
errors = []
dofs = []

for alpha in alphas:
    for theta in thetas:
        # run fixed-theta AMR
        action = np.array([theta])
        actions.append(action.item())
        print("action = ", action.item())
        env.reset(new_alpha = False)
        env.alpha = alpha
        done = False
        episode_cost_tmp = 0
        errors_tmp = [env.global_error]
        dofs_tmp = [env.sum_of_dofs]
        steps_taken = 0
        
        while not done:
            _, reward, done, info = env.step(action)
            if steps_taken > num_iterations:
                print("*** fixed action exceeded max step threshold of ", num_iterations, "steps.")
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

        # find final cumulative degrees of freedom
        cum_dofs = np.cumsum(dofs[-1])[-1]

        # print results to file
        result_str = str(alpha) + ", " + str(theta) + ", " + str(cum_dofs) + ", " + str(errors[-1][-1]) + "\n"
        output_file.write(result_str)

# close output file
output_file.close()
