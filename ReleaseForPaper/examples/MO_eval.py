import numpy as np

def evaluate_alpha(env, trainer, alpha):
    env.alpha = alpha
    print("evaluating alpha = {}".format(alpha))

    ## Enact trained policy 
    env.trainingmode = False
    rlactions = []
    obs = env.reset(new_alpha = False)

    done = False
    rlepisode_cost = 0
    rlerrors = [env.global_error]
    rldofs = [env.sum_of_dofs]
    
    while not done:
        action = trainer.compute_single_action(obs, explore=False)
        obs, reward, done, info = env.step(action)

        rlactions.append(action[0])
        rlepisode_cost -= reward
        print("step = ", env.k)
        print("action = ", action.item())
        print("Num. Elems. = ", env.mesh.GetNE())
        print("episode cost = ", rlepisode_cost)
        rldofs.append(info['num_dofs'])
        rlerrors.append(info['global_error'])

    # save final errors in file for each alpha
    cum_dofs = np.cumsum(rldofs)[-1]
    error = rlerrors[-1]
    
    return cum_dofs, error

def MO_eval(env, trainer):
    '''
    
    '''

    ##  evaluate for alpha in {0.1, 0.2, ..., 0.9}
    alphas_to_eval = 0.1*np.array(range(1, 10, 1))
    cum_dofs       = []
    errors         = []

    print("starting multi-objective evaluation")
    for alpha in alphas_to_eval:
        alpha_cum_dofs, alpha_error = evaluate_alpha(env, trainer, alpha)

        # save final errors in file for each alpha
        cum_dofs.append(alpha_cum_dofs)
        errors.append(alpha_error)

    env.trainingmode = True # set training mode back for next batch
    
    ## sort data by increasing cumulative dofs
    sorted_indices = np.argsort(cum_dofs)
    cum_dofs       = np.array(cum_dofs)[sorted_indices]
    errors         = np.array(errors)[sorted_indices] 

    ## approximate area under cumulative dofs vs error curve using Riemann sum of log of data
    print("cum_dofs = {}".format(cum_dofs))
    widths = cum_dofs[0] * np.ones(len(cum_dofs))
    for i in range(0, len(cum_dofs) - 1):
        widths[i+1] = np.log(cum_dofs[i+1]) - np.log(cum_dofs[i])
    heights = errors # we don't take the log of errors because we don't want negative numbers
    rect_areas = np.multiply(widths, heights) # multiply widths and heights element-wise

    # add up rect_areas to approximate area under dofs vs error curve
    area = np.sum(rect_areas)

    ## return area as a metric for how well training is going
    return area
