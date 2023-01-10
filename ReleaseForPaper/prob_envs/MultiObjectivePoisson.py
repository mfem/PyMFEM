from prob_envs.Poisson import Poisson
import numpy as np
from utils.Statistics import Statistics, GlobalError
from gym import spaces
from utils.ReentrantCorner import *
from mfem.ser import intArray
from os.path import expanduser, join
from utils.Solution_LShaped import *
from utils.Solution_Wavefront import *
from utils.Solution_SinSin import *
import os
import ray

class MultiObjPoisson(Poisson):
    '''
    This class inherits from the Poisson class, but allows for cost functions which are a linear combination of the 
    cost due to DOFs and the cost due to global error.
    '''

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.optimization_type  = kwargs.get('optimization_type','multi_objective')
        self.alpha              = kwargs.get('alpha', 0.5)              # default is to weight each objective equally
        self.observe_alpha      = kwargs.get('observe_alpha', True)     # observe alpha, so policy will depend on alpha
        self.observe_budget     = kwargs.get('observe_budget', True)    # observe budget = (current step #)/(total number of steps)
        self.num_iterations     = kwargs.get('num_iterations', 10)      # decide what a good default number of iterations is

        # define range of expected values for the observations
        if self.observe_alpha == True and self.observe_budget == True:
            self.observation_space = spaces.Box(low = np.array([-np.inf,-np.inf, 0.0, 0.0]), high= np.array([np.inf, np.inf, 1.0, 1.0]))
        elif self.observe_alpha == False and self.observe_budget == False:
            self.observation_space = spaces.Box(low = np.array([-np.inf,-np.inf]), high= np.array([np.inf, np.inf]))
        else: 
            # this is the case where only one of observe_alpha and observe_budget are True
            self.observation_space = spaces.Box(low = np.array([-np.inf,-np.inf, 0.0]), high= np.array([np.inf, np.inf, 1.0]))

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
            elif self.sum_of_dofs >= self.dof_threshold:
                done = True
                print("*** dof threshold was exceeded ***")
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

            # define observation, which depends on observe_alpha and observe_budget
            obs = [stats.mean, stats.variance]
            if self.observe_alpha == True: 
                obs.append(self.alpha)

            if self.observe_budget == True:
                budget = self.k/self.num_iterations
                obs.append(budget)

            return np.array(obs)
        else:
            return super().GetObservation(self)

    def reset(self, new_alpha = True):
        # choose new random alpha if alpha is included as an observation
        if self.optimization_type == 'multi_objective':
            if self.observe_alpha == True and new_alpha == True: 
                # set random alpha value
                self.alpha = np.random.uniform(low = 0., high = 1.)
        
        return super().reset()

class ADF_MultiObjPoisson(Poisson):
    '''
    This class inherits from the Poisson class, but allows for cost functions which are given by
    cost = (cost due to cumulative dofs) * |tau - (cost due to final error)|/delta, 
    where tau is a target error and delta is a positive parameter that defines a sort of "width".
    '''

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.optimization_type  = kwargs.get('optimization_type','multi_objective')
        self.observe_budget     = kwargs.get('observe_budget', True)    # observe budget = (current step #)/(total number of steps)
        self.num_iterations     = kwargs.get('num_iterations', 10)      # decide what a good default number of iterations is
        self.tau_min            = kwargs.get('tau_min', np.log2(10**-4))# minimium target error to train for
        self.N_anneal           = kwargs.get('N_anneal', 50)            # number of target errors between tau_min and tau_max to train for
        self.M_warm             = kwargs.get('M_warm',   50)            # number of batches in warming phase
        self.M_anneal           = kwargs.get('M_anneal', 20)            # number of batches for each tau in annealing phase
        self.batch_size         = kwargs.get('batch_size', 100)         # number of episodes per batch
        
        # define range of expected values for the observations
        if self.observe_budget == True:
            # observation space includes mean, variance, tau, budget
            self.observation_space = spaces.Box(low = np.array([-np.inf,-np.inf, -np.inf, 0.0]), high= np.array([np.inf, np.inf, np.inf, 1.0]))
        else: 
            # observation space includes mean, variance, tau
            self.observation_space = spaces.Box(low = np.array([-np.inf,-np.inf, -np.inf]), high= np.array([np.inf, np.inf, np.inf]))

        self.ADF_Params = ray.get_actor("parameters")

        # arbitrarily set tau_max to start (this should be reset once we call FindParameters)
        self.reset_tau = False
        self.tau_max = 0.1
        
    def FindParameters(self):
        self.mesh = mfem.Mesh(self.initial_mesh)
        self.Setup()
        self.AssembleAndSolve()
        self.errors = self.GetLocalErrors()
        self.global_error = max(GlobalError(self.errors),1e-12)
        error_init = np.log2(self.global_error)

        tau_max = np.log2(0.9) + error_init; # set tau_max to log(0.9 * initial error) using log rule
        # tau_max      = self.tau_min + self.N_anneal/(self.N_anneal + 1)*(error_init - self.tau_min)
        if self.N_anneal != 0:
            tau_step  = (tau_max - self.tau_min)/self.N_anneal
        else:
            tau_step = 0
        delta_warm   = (tau_max - self.tau_min)/2
        delta_anneal = (tau_max - self.tau_min)/10

        self.tau_max = tau_max; 
        #self.reset_tau = True;
        #print('reset tau set to True')
        return tau_max, tau_step, delta_warm, delta_anneal

    def step(self, action):
        self.reset_tau = False  # maybe want this to be True if adjusting Tau **************
        if self.optimization_type == 'multi_objective':
            self.k += 1 # increment the step index
            self.UpdateMesh(action)
            #print('Theta = {}'.format(action))
            # find errors and num dofs
            self.AssembleAndSolve()
            self.errors = self.GetLocalErrors()
            num_dofs = self.fespace.GetTrueVSize()
            global_error = GlobalError(self.errors)

            # update cost = alpha*(dof cost) + (1-alpha)*(error cost)
            if True:
                dofs_cost = np.log2(self.sum_of_dofs + num_dofs);
                error_cost = np.log2(global_error);
            else: 
                dofs_cost  = np.log2(1.0 + num_dofs/self.sum_of_dofs)
                error_cost = np.log2(global_error/self.global_error)

           # cost = dofs_cost * np.abs(self.tau - error_cost)/self.delta

            self.sum_of_dofs += num_dofs
            self.global_error = global_error

            if self.k >= self.num_iterations or self.sum_of_dofs >= self.dof_threshold:
                ## the cost computation in the next line penalizes the difference of the logs of the error:
                # cost = dofs_cost * np.abs(ray.get(self.ADF_Params.get_tau.remote()) - error_cost)/ray.get(self.ADF_Params.get_delta.remote())

                ## the cost computation in next line penalizes the log of the difference of the error
                cost = dofs_cost * np.abs(ray.get(self.ADF_Params.get_tau.remote()) - error_cost)/ray.get(self.ADF_Params.get_delta.remote())
                done = True
                # if self.sum_of_dofs >= self.dof_threshold:
                #     print("*** dof threshold was exceeded: penalizing cost ***")
                #     cost = 100*cost
                #     print("cost is ", cost)
                # tau_used = ray.get(self.ADF_Params.get_tau.remote())
                # print("tau used = {}".format(tau_used))
                # print("dofs cost = {}, F_2 = {}".format(dofs_cost, np.abs(tau_used - error_cost)/self.delta))
                # print("sum of dofs at stop  = {}".format(self.sum_of_dofs), flush=True)
                # print("global error at stop = {}".format(self.global_error), flush=True)
                print("at stop: error=", np.round(self.global_error,4), "cumdofs=", self.sum_of_dofs, "cost=", cost, 
                        "tau diff=", np.round(ray.get(self.ADF_Params.get_tau.remote()) - error_cost,4), "tau=", np.round(ray.get(self.ADF_Params.get_tau.remote()), 4) ,"delta=", np.round(ray.get(self.ADF_Params.get_delta.remote()),4))
            else:
                cost = 0
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
            stats = Statistics(self.errors, num_dofs = num_dofs)

            # define observation, which depends on observe_alpha and observe_budget
            obs = [stats.mean, stats.variance, ray.get(self.ADF_Params.get_tau.remote())]

            if self.observe_budget == True:
                budget = self.k/self.num_iterations
                obs.append(budget)

            return np.array(obs)
        else:
            return super().GetObservation(self)

    def reset(self):
        if self.reset_tau:
        # reset tau to random value between tau_min and tau_max
            tau = np.random.uniform(low = self.tau_min, high = ray.get(self.ADF_Params.get_tau_init.remote()))
            #print("Tau_min = {}, Tau_max = {}.".format(self.tau_min, ray.get(self.ADF_Params.get_tau_init.remote())))
            self.ADF_Params.set_tau.remote(tau) 
            #print("Reset tau to {}".format(tau))
        return super().reset()

class Angle_MultiObjPoisson(MultiObjPoisson):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)  
        problem_type = kwargs.get('problem_type','lshaped')
        if (problem_type == 'lshaped'):
            delattr(self, 'BC')
            self.BC = mfem.NumbaFunction(ReentrantCornerExact, 2, True).GenerateCoefficient()
        elif (problem_type == 'wavefront'):
            print("*** Should not be using RobustAnglePoisson class with wavefront problem - exiting ***")
            exit()
        else:
            delattr(self, 'BC')
            self.BC = mfem.ConstantCoefficient(0.0)
        self.num_unif_ref = kwargs.get('num_unif_ref',1)
        self.Lshapedmeshfile = expanduser(join(os.path.dirname(__file__), '../..', 'data', 'l-shape-benchmark.mesh'))
        self.circlemeshfile = expanduser(join(os.path.dirname(__file__), '../..', 'data', 'circle_3_4.mesh'))

        self.angle = kwargs.get('angle_lower', np.pi * 0.5)
        self.angle_lower = kwargs.get('angle_lower', np.pi * 0.25)
        self.angle_upper = kwargs.get('angle_upper', np.pi * 0.75)
        if (problem_type == 'lshaped'):
            self.set_angle(self.angle)

    def set_angle(self, angle):
        # print("Setting env angle to ", self.angle, flush=True)
        self.angle = angle
        self.BC.SetTime(self.angle)

        if self.mesh_name =='l-shape-benchmark.mesh':
            self.initial_mesh = ReentrantCornerMesh(self.angle, self.Lshapedmeshfile)
        elif self.mesh_name =='circle_3_4.mesh':
            self.initial_mesh = ReentrantCornerMesh(self.angle, self.circlemeshfile)
        else:
            print("Something went wrong - see set_angle routine; exiting")
            exit()
        for _ in range(self.num_unif_ref):
            self.initial_mesh.UniformRefinement()
        self.initial_mesh.EnsureNCMesh()

    def reset(self, random_angle=True, new_alpha = True):

        if random_angle:
            angle = np.random.uniform(self.angle_lower, self.angle_upper, 1).item()
            # print("Resetting env angle to ", angle)
            self.set_angle(angle)
        return super().reset(new_alpha = new_alpha)




class hp_Angle_MultiObjPoisson(Angle_MultiObjPoisson):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)  
        self.action_space = spaces.Box(low = np.array([0.0, 0.0]), high = np.array([0.999, 0.999]), shape=(2,), dtype = np.float32)

    def Prefine(self, theta, rho):   
        mark_to_p_refine = []
        threshold = theta * np.max(self.errors)
        for i in range(self.mesh.GetNE()):
            if threshold >= self.errors[i]:
                mark_to_p_refine.append((i, self.errors[i]))
        mark_to_p_refine.sort(key=lambda x:x[1], reverse=True)
        for i in range(len(mark_to_p_refine)):
            if mark_to_p_refine[i][1] > rho * np.max(self.errors):
                current_element = mark_to_p_refine[i][0]
                current_order = self.fespace.GetElementOrder(current_element)
                self.fespace.SetElementOrder(current_element, current_order + 1)
        
        self.fespace.Update(False)
        self.x.Update()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        # self.fespace.UpdatesFinished()
        self.a.Update()
        self.b.Update()

    def CloseMesh(self, delta_p = 1):
        # Loop through all elements in mesh until the maximum difference in polynomial
        # orders across all edges is no more than delta_p
        neighbor_table = self.mesh.ElementToElementTable()
        while True:
            mesh_closed = True
            elements_to_p_refine = []
            for i in range(self.mesh.GetNE()):
                neighbor_row = neighbor_table.GetRow(i)
                row_size = neighbor_table.RowSize(i)
                neighbor_array = intArray(row_size)
                neighbor_array.Assign(neighbor_row)
                for l in range(row_size):
                    neighbor_order = self.fespace.GetElementOrder(neighbor_array[l])
                    if neighbor_order - self.fespace.GetElementOrder(i) > delta_p:
                        elements_to_p_refine.append(i)
                        mesh_closed = False
            p_refine_elements = np.unique(elements_to_p_refine).tolist()
            for k in range(len(p_refine_elements)):
                current_element = p_refine_elements[k]
                current_order = self.fespace.GetElementOrder(current_element)
                self.fespace.SetElementOrder(current_element, current_order + 1)
            if mesh_closed:
                break

        self.fespace.Update(False)
        self.x.Update()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        # self.fespace.UpdatesFinished()
        self.a.Update()
        self.b.Update()

    # overriding UpdateMesh in Poisson (grand)-parent class
    def UpdateMesh(self, action):
        action = np.clip(action, 0.0, 1.0)
        theta = action[0].item() 
        rho = action[1] * theta 
        # self.refinement_strategy == 'max'
        self.Prefine(theta, rho)
        self.Refine(theta)
        self.CloseMesh()
    
    def GetElementVertices(self, k):
        Tr = self.mesh.GetElementTransformation(k)
        physical_pts = np.zeros((4,2))
        reference_pt = mfem.IntegrationPoint()
        for i in range(2):
            for j in range(2):
                reference_pt.Set(float(i),float(j),0.0,0.0)
                physical_pts[i+2*j,:] = Tr.Transform(reference_pt)
        return physical_pts

    def RenderHPmesh(self, gfname=None):
        ordersfec = mfem.L2_FECollection(0, self.dim)
        ordersfes = mfem.FiniteElementSpace(self.mesh, ordersfec)
        orders = mfem.GridFunction(ordersfes)
        for i in range(0, self.mesh.GetNE()):
            elem_dofs = 0
            elem_dofs = ordersfes.GetElementDofs(i)
            # orders[elem_dofs[0]] = self.errors[i]
            orders[elem_dofs[0]] = self.fespace.GetElementOrder(i)
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock.send_solution(self.mesh, orders)
        title = "step " + str(self.k)
        sol_sock.send_text('keys ARjlmp*******' + " window_title '" + title)
        sol_sock.send_text("valuerange 1.0 8.0 \n")
        if gfname:
            orders.Save(str(gfname), 8)  # second input is "precision"
