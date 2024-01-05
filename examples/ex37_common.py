'''
   PyMFEM example 37 - Serial/Parallel Shared code
'''

import os
from os.path import expanduser, join

from numpy import sqrt, log, exp
import numpy as np
from numba import njit
from numba.types import float64

from mfem import mfem_mode

if mfem_mode == 'serial':
   import mfem.ser as mfem
   from mfem.ser import intArray, doubleArray
   use_parallel = False
else:
   import mfem.par as mfem
   from mfem.par import intArray, doubleArray
   use_parallel = True   

@njit(float64(float64))
def inv_sigmoid(x):
    '''
    Inverse sigmoid function
    '''
    tol = 1e-12;
    x = min(max(tol,x),1.0-tol)
    return log(x/(1.0-x))


@njit(float64(float64))
def sigmoid(x):
    '''
    Sigmoid function
    '''
    if x >= 0:
        return 1.0/(1.0 + exp(-x))
    else:
       return exp(x)/(1.0 + exp(x))

@njit(float64(float64))    
def der_sigmoid(x):
    '''
    Derivative of sigmoid function    
    '''
    tmp = sigmoid(-x)
    return tmp - tmp**2


def MappedGridFunctionCoefficient(gf, func):

    c_gf = mfem.GridFunctionCoefficient(gf)
    @mfem.jit.scalar(dependency=(c_gf,))
    def coeff(ptx, c_gf):
        return func(c_gf)
    return coeff

def DiffMappedGridFunctionCoefficient(gf, other_gf, func, comp=1):
    c_gf = mfem.GridFunctionCoefficient(gf)
    c_ogf = mfem.GridFunctionCoefficient(other_gf)
    
    @mfem.jit.scalar(dependency=(c_gf, c_ogf))
    def coeff(ptx, c_gf, c_ogf):
        return func(c_gf) - func(c_ogf)
    return coeff    
    
class SIMPInterpolationCoefficient():
    '''
    Python Note: 
       Approach here is to replace Eval in C++ example using the dependency 
       feature of mfem.jit.

       In order to avoid repeating Numba-Jitting in iteration loop, we use
       SetGridFunction to update the GridFunction referred from 
       GridFunctionCoefficient.
    '''
    def __init__(self, rho_filter, min_val=1e-6, max_val = 1.0, exponent=3):
       val = mfem.GridFunctionCoefficient(rho_filter)
       @mfem.jit.scalar(dependency=(val,))
       def coeff(ptx, val):
           coeff = min_val + val**exponent*(max_val-min_val)
           return coeff

       self.c_gf = val
       self.coeff = coeff

    def Update(self, rho_filter):
       self.c_gf.SetGridFunction(rho_filter)

class StrainEnergyDensityCoefficient():
    '''
    Python Note: 
       Approach here is to replace Eval and GetVectorGradient method call in C++
       using the dependency feature of mfem.jit.

       GetVectorGradient is mimiced by creating GradientGridFunctionCoefficient
       for each component of u vector. Note GridFunction(fes, u.GetDataArray())
       reuses the data array from u.
    '''
    def __init__(self, llambda, mu, u, rho_filter, rho_min=1e-6, exponent=3.0):    
       assert rho_min >= 0.0, "rho_min must be >= 0"
       assert rho_min < 1.0,  "rho_min must be > 1"

       fes = u.FESpace()
       assert fes.GetOrdering() == mfem.Ordering.byNODES, "u has to use byNODES ordering"

       mesh = fes.GetMesh()
       dim = mesh.Dimension()
       assert dim == 2, "dim must be two."

       fec = fes.FEColl()
       fes = mfem.FiniteElementSpace(mesh, fec)
       size = len(u.GetDataArray())

       u1 = mfem.GridFunction(fes, mfem.Vector(u.GetDataArray()))   # first component 
       u2 = mfem.GridFunction(fes, mfem.Vector(u.GetDataArray()), size//2) # second component

       c_gradu1 = mfem.GradientGridFunctionCoefficient(u1)
       c_gradu2 = mfem.GradientGridFunctionCoefficient(u2)


       c_rho_filter = mfem.GridFunctionCoefficient(rho_filter)
       @mfem.jit.scalar(dependency=(llambda, mu, c_gradu1, c_gradu2, c_rho_filter))
       def coeff(ptx, L, M, grad1, grad2, val):
           div_u = grad1[0] + grad2[1]
           density = L*div_u*div_u

           grad = np.zeros(shape=(2,2), dtype=np.float64)
           grad[0,0] = grad1[0]
           grad[0,1] = grad1[1]        
           grad[1,0] = grad2[0]
           grad[1,1] = grad2[1]

           for i in range(2):
               for j in range(2):
                  density += M*grad[i, j]*(grad[i, j]+ grad[j,i])
           return -exponent * val**(exponent-1.0)*(1-rho_min)*density

       self.fes = fes
       self.size = size
       self.u1u2 = (u1, u2)
       self.dependency = (c_gradu1, c_gradu2, c_rho_filter)
       self.coeff = coeff
        
    def Update(self, u, rho_filter):
       u1 = mfem.GridFunction(self.fes, mfem.Vector(u.GetDataArray()))
       u2 = mfem.GridFunction(self.fes, mfem.Vector(u.GetDataArray()),
                              self.size//2)
       self.dependency[0].SetGridFunction(u1)
       self.dependency[1].SetGridFunction(u2)
       self.dependency[2].SetGridFunction(rho_filter)
       self.u1u2 = (u1, u2)       
    
def VolumeForceCoefficient(r, center, force):

    @mfem.jit.vector(shape=(len(center),))
    def coeff(ptx):
        cr = sqrt(sum((ptx - center)**2))
        if cr < r:
            return np.array((force[0], force[1]))
        else:
            return np.array((0.0, 0.0))
    return coeff

class DiffusionSolver():
    def __init__(self):
        self.rhscf = None
        self.neumann_cf = None
        self.masscf = None
        self.essbdr_cf = None
        self.gradient_cf = None        

    def SetupFEM(self):
        dim = self.mesh.Dimension()
        self.fec = mfem.H1_FECollection(self.order, dim)
        self.fes = mfem.FiniteElementSpace(self.mesh, self.fec)

        if self.ess_bdr.Size() == 0 and mesh.bdr_attributes.Size() > 0:
           self.ess_bdr = mfem.intArray([1]*mesh.bdr_attributes.Max())
            
    def Solve(self):
        A = mfem.OperatorPtr()
        B = mfem.Vector()
        X = mfem.Vector()
        ess_tdof_list = mfem.intArray()
        self.fes.GetEssentialTrueDofs(self.ess_bdr, ess_tdof_list)

        self.u = mfem.GridFunction(self.fes)
        self.u.Assign(0.0)

        b = mfem.LinearForm(self.fes)

        if self.rhscf is not None:
           itg = mfem.DomainLFIntegrator(self.rhscf)
           b.AddDomainIntegrator(itg)

        if self.neumann_cf is not None:
            assert self.neumann_bdr.Size()>0, "neumann_bdr attributes not provided"
            b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(self.neumann_cf),
                                    self.neumann_bdr)
        elif self.gradient_cf is not None:
            assert self.neumann_bdr.Size()>0, "neumann_bdr attributes not provided"
            b.AddBoundaryIntegrator(mfem.BoundaryNormalLFIntegrator(self.gradient_cf),
                                    self.neumann_bdr)

        b.Assemble()

        a = mfem.BilinearForm(self.fes)
        a.AddDomainIntegrator(mfem.DiffusionIntegrator(self.diffcf))
        if self.masscf is not None:
            a.AddDomainIntegrator(mfem.MassIntegrator(self.masscf))
        a.Assemble()
        if self.essbdr_cf is not None:
            self.u.ProjectBdrCoefficient(essbdr_cf,ess_bdr)

        a.FormLinearSystem(ess_tdof_list, self.u, b, A, X, B)
        AA = A.AsSparseMatrix()
        M = mfem.GSSmoother(AA)
        cg = mfem.CGSolver()
        cg.SetRelTol(1e-12)
        cg.SetMaxIter(10000)
        cg.SetPrintLevel(0)
        cg.SetPreconditioner(M)
        cg.SetOperator(A)
        cg.Mult(B, X)
        
        a.RecoverFEMSolution(X, b, self.u)
        self.b = b

    def GetFEMSolution(self):
        return self.u
    def GetLinearForm(self):
        return self.b
    def SetRHSCoefficient(self, rhscf):
        self.rhscf = rhscf
        
class LinearElasticitySolver():
    def __init__(self):
        self.rhs_cf = None
        self.essbdr_cf = None


    def SetupFEM(self):
        dim = self.mesh.Dimension()
        self.fec = mfem.H1_FECollection(self.order, dim, mfem.BasisType.Positive)
        self.fes = mfem.FiniteElementSpace(self.mesh, self.fec, dim)
        self.u = mfem.GridFunction(self.fes)
        self.u.Assign(0.0)

    def Solve(self):
        A = mfem.OperatorPtr()
        B = mfem.Vector()
        X = mfem.Vector()
        
        ess_tdof_list = mfem.intArray()
        self.fes.GetEssentialTrueDofs(self.ess_bdr,ess_tdof_list)
        
        x = mfem.GridFunction(self.fes)
        x .Assign(0.0)
        
        self.u.Assign(0.0)
        b = mfem.LinearForm(self.fes)

        if self.rhs_cf is not None:
            b.AddDomainIntegrator(mfem.VectorDomainLFIntegrator(self.rhs_cf))
        b.Assemble()

        a = mfem.BilinearForm(self.fes)
        a.AddDomainIntegrator(mfem.ElasticityIntegrator(self.lambda_cf, self.mu_cf))
        a.Assemble()
        if self.essbdr_cf is not None:
             u.ProjectBdrCoefficient(self.essbdr_cf, self.ess_bdr)

        a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)

        AA = A.AsSparseMatrix()                                     
        M = mfem.GSSmoother(AA)
        cg = mfem.CGSolver()

        cg.SetRelTol(1e-10)
        cg.SetMaxIter(10000)
        cg.SetPrintLevel(0)
        cg.SetPreconditioner(M)
        cg.SetOperator(A)
        cg.Mult(B, X)
        a.RecoverFEMSolution(X, b, x)

        self.u += x
        self.b = b

    def GetFEMSolution(self):
        return self.u
    def GetLinearForm(self):
        return self.b
        

