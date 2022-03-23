from mfem._ser.coefficient import ConstantCoefficient
import os
from os.path import expanduser, join
import gym
from gym import spaces
import numpy as np
from tensorflow.python.ops.gen_array_ops import Const
import mfem.ser as mfem
from mfem.ser import intArray
from utils.Statistics import Statistics, GlobalError
from utils.Solution_LShaped import *

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import dsolve
def SolveSparseSystem(A,b,x):
   A_I = A.GetIArray()
   A_J = A.GetJArray()
   data = A.GetDataArray()
   n = A.Height()
   SpA = csr_matrix( (data, A_J, A_I), shape=(n, n) )
   npX = np.empty([n])
   for j in range(n):
      npX[j] = b[j]
   Y = dsolve.spsolve(SpA, npX, use_umfpack=False)
   for j in range(n):
      x[j] = Y[j]


class Poisson(gym.Env):

    def __init__(self,**kwargs):
        super().__init__()
        self.nc_limit = 1 # maximum number of hanging nodes per element = 2^nc_limit - 1 
        self.order = kwargs.get('order',2)
        
        problem_type = kwargs.get('problem_type','lshaped')
        if (problem_type == 'lshaped'):
            self.BC = mfem.NumbaFunction(LShapedExact, 2).GenerateCoefficient()
            self.RHS = mfem.ConstantCoefficient(0.0)
            self.coeff = mfem.ConstantCoefficient(1.0)
        else:
            self.BC = mfem.ConstantCoefficient(0.0)
            self.RHS = mfem.ConstantCoefficient(1.0)
            self.coeff = mfem.ConstantCoefficient(1.0)
        
        self.optimization_type = kwargs.get('optimization_type','error_threshold')
        self.error_threshold = kwargs.get('error_threshold',1e-3)
        self.dof_threshold = kwargs.get('dof_threshold',5e4)

        mesh_name = kwargs.get('mesh_name','l-shape-benchmark.mesh')
        num_unif_ref = kwargs.get('num_unif_ref',1)
        self.meshfile = expanduser(join(os.path.dirname(__file__), '../..', 'data', mesh_name))
        mesh = mfem.Mesh(self.meshfile)
        mesh.EnsureNCMesh()
        for _ in range(num_unif_ref):
            mesh.UniformRefinement()
        self.dim = mesh.Dimension()
        self.initial_mesh = mesh
        
        self.action_space = spaces.Box(low=0.0, high=0.999, shape=(1,), dtype=np.float32)
        self.observation_space = spaces.Box(low = np.array([0.0,-np.inf,-np.inf]), 
                                            high= np.array([1.0, np.inf, np.inf]))
        # self.observation_space = spaces.Box(low=-np.inf, high=np.inf, shape=(3,))
        # self.observation_space = spaces.Box(low=-np.inf, high=np.inf, shape=(5,))
        self.error_file_name = kwargs.get('error_file_name','./ReleaseForPaper/out/errors.csv')

        self.zero = mfem.ConstantCoefficient(0.0)
        self.zero_vector = mfem.Vector(self.dim)
        self.zero_vector.Assign(0.0)
        self.zerovector = mfem.VectorConstantCoefficient(self.zero_vector)

        self.global_errors = [ [] for _ in range(4)]
        self.trainingmode = True

    
    def reset(self):
        self.k = 0
        self.mesh = mfem.Mesh(self.initial_mesh)
        self.Setup()
        self.AssembleAndSolve()
        self.errors = self.GetLocalErrors()
        self.global_error = max(GlobalError(self.errors),1e-12)
        self.sum_of_dofs = self.fespace.GetTrueVSize()
        obs = self.GetObservation()
        return obs
    
    def step(self, action):
        self.k += 1
        self.UpdateMesh(action)
        if self.optimization_type == 'error_threshold':
            self.AssembleAndSolve()
            self.errors = self.GetLocalErrors()
            num_dofs = self.fespace.GetTrueVSize()
            self.global_error = GlobalError(self.errors)
            cost = np.log2(1.0 + num_dofs/self.sum_of_dofs)
            self.sum_of_dofs += num_dofs
            if self.global_error < self.error_threshold:
                done = True
            else:
                done = False
            if self.trainingmode:
                if self.sum_of_dofs > self.dof_threshold:
                    cost += 10.0
                    done = True
        elif self.optimization_type == 'dof_threshold':
            num_dofs = self.fespace.GetTrueVSize()
            self.sum_of_dofs += num_dofs
            if self.sum_of_dofs > self.dof_threshold:
                cost = 0.0
                done = True
            else:
                self.AssembleAndSolve()
                self.errors = self.GetLocalErrors()
                global_error = GlobalError(self.errors)
                cost = np.log2(global_error/self.global_error)
                self.global_error = global_error
                done = False
        else:
            assert False, "Houston we've got a problem"
        if done == False:
            obs = self.GetObservation()
        else:
            obs = np.zeros_like(self.GetObservation())
        info = {'global_error':self.global_error, 'num_dofs':num_dofs, 'max_local_errors':np.amax(self.errors)}
        return obs, -cost, done, info
    
    def render(self):
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock.send_solution(self.mesh,  self.x)
        title = "step " + str(self.k)
        sol_sock.send_text("window_title '" + title)

    def Setup(self):
        dim = self.mesh.Dimension()
        fec = mfem.H1_FECollection(self.order, dim)
        self.fespace = mfem.FiniteElementSpace(self.mesh, fec)
        self.a = mfem.BilinearForm(self.fespace)
        self.b = mfem.LinearForm(self.fespace)
        integ = mfem.DiffusionIntegrator(self.coeff)
        self.a.AddDomainIntegrator(integ)
        self.b.AddDomainIntegrator(mfem.DomainLFIntegrator(self.RHS))
        self.x = mfem.GridFunction(self.fespace)
        self.x.Assign(0.0)
        self.ess_bdr = intArray(self.mesh.bdr_attributes.Max())
        self.ess_bdr.Assign(1)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        self.flux_fespace = mfem.FiniteElementSpace(self.mesh, fec, dim)
        self.estimator =  mfem.NewZienkiewiczZhuEstimator(integ, self.x, self.flux_fespace,
                                                          own_flux_fes = False)

    def GetObservation(self):
        num_dofs = self.fespace.GetTrueVSize()
        stats = Statistics(self.errors, num_dofs=num_dofs)
        if self.optimization_type == 'error_threshold':
            # budget = -np.log( abs(self.error_threshold - self.global_error)/self.global_error + 1e-16)
            budget = self.error_threshold/self.global_error
        else:
            # budget = np.log( self.sum_of_dofs / (abs(self.dof_threshold - self.sum_of_dofs) + 1e-12 ))
            budget = self.sum_of_dofs/self.dof_threshold
        
        obs = [budget, stats.mean, stats.variance]
        return np.array(obs)

    def AssembleAndSolve(self):
        self.a.Assemble()
        self.b.Assemble()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        ess_tdof_list = intArray()
        self.fespace.GetEssentialTrueDofs(self.ess_bdr, ess_tdof_list)
        A = mfem.OperatorPtr()
        B = mfem.Vector();  X = mfem.Vector()
        self.a.FormLinearSystem(ess_tdof_list, self.x, self.b, A, X, B, 1)
        AA = mfem.OperatorHandle2SparseMatrix(A)
        # M = mfem.GSSmoother(AA)
        # mfem.PCG(AA, M, B, X, -1, 2000, 1e-12, 0.0)
        AA.Threshold(0.0, False)
        SolveSparseSystem(AA,B,X)
        self.a.RecoverFEMSolution(X,self.b,self.x)
        self.solution_norm = self.x.ComputeGradError(self.zerovector) + 1e-12

    def GetLocalErrors(self):
        self.estimator.Reset()
        self.mfem_errors = self.estimator.GetLocalErrors()
        errors = np.zeros(self.mesh.GetNE())
        for i in range(self.mesh.GetNE()):
            self.mfem_errors[i] /= self.solution_norm
            errors[i] = self.mfem_errors[i]
        return errors

    def RenderMesh(self):
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        zerogf = mfem.GridFunction(self.fespace)
        zerogf.Assign(0.0)
        sol_sock.send_solution(self.mesh, zerogf)
        title = "step " + str(self.k)
        sol_sock.send_text('keys ARjlmp*******' + " window_title '" + title)

    def UpdateMesh(self, action):
        action = np.clip(action, 0.0, 1.0)
        theta = action.item() # refinement threshold
        self.Refine(theta)

    # "Refine" does h-refine only:
    #   input theta is in [0,1]
    #   h-refine all elements with error geq theta * (max error)
    def Refine(self, theta):
        threshold = theta * np.max(self.errors)
        self.mesh.RefineByError(self.mfem_errors,threshold, -1, self.nc_limit)
        self.fespace.Update()
        self.x.Update()
        self.a.Update()
        self.b.Update()