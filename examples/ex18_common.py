class InitialCondition(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
            
class FE_Evolution(mfem.TimeDependentOperator):
    pass       
class DomainIntegratoer(mfem.BilinearFormIntegrator):
    pass                 
class FaceIntegratoer(mfem.NonlinearFormIntegrator):
    pass
class RiemanSovler():
    pass
def StateIsPhysical(sdata, dim):
    pass
def ComputePressure(state, dim):
    pass
def ComputeFlux(state, dim):
    pass
def ComputeFluxDotN(state, nor, fluxN):
    pass
def ComputeMaxCharSpeed(state, dim):
    pass
