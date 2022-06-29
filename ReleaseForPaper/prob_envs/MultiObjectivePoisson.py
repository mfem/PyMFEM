from prob_envs.Poisson import Poisson
import numpy as np
from utils.Statistics import Statistics, GlobalError
from gym import spaces

class MultiObjPoisson(Poisson):
    '''
    This class inherits from the Poisson class, but allows for cost functions which are a linear combination of the 
    cost due to DOFs and the cost due to global error.
    '''

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.optimization_type = kwargs.get('optimization_type','multi_objective')
        self.alpha = kwargs.get('alpha', 0.5) # default is to weight each objective equally
        self.num_iterations = kwargs.get('num_iterations', 10) # decide what a good default number of iterations is
        self.observation_space = spaces.Box(low = np.array([-np.inf,-np.inf]), high= np.array([np.inf, np.inf]))

    def step(self, action):
        if self.optimization_type == 'multi_objective':
            self.k += 1 # increment the step index
            self.UpdateMesh(action)

            # find errors and num dofs
            self.AssembleAndSolve()
            self.errors = self.GetLocalErrors()
            num_dofs = self.fespace.GetTrueVSize()
            global_error = GlobalError(self.errors)

            # update cost = alpha*(dof cost) + (1-alpha)*(error cost)
            if self.k == 1:
                cost = self.alpha * np.log2(self.sum_of_dofs + num_dofs) + (1 - self.alpha) * np.log2(global_error)
            else: 
                cost = self.alpha * np.log2(1.0 + num_dofs/self.sum_of_dofs) + (1 - self.alpha) * np.log2(global_error/self.global_error)

            self.sum_of_dofs += num_dofs
            self.global_error = global_error

            if self.k >= self.num_iterations:
                done = True
            else:
                done = False

            if done == False:
                obs = self.GetObservation()
            else:
                obs = np.zeros_like(self.GetObservation())

            info = {'global_error':self.global_error, 'num_dofs':num_dofs, 'max_local_errors':np.amax(self.errors)}
            return obs, -cost, done, info

        # if using a single objective, use the parent class
        else: 
            return super().step(self, action)

    def GetObservation(self):
        if self.optimization_type == 'multi_objective':
            num_dofs = self.fespace.GetTrueVSize()
            stats = Statistics(self.errors, num_dofs=num_dofs)
            obs = [stats.mean, stats.variance]
            return np.array(obs)

        else:
            return super().GetObservation(self)
