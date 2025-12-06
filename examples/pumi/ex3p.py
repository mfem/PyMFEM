'''
   MFEM example 3 parallel version / using pumi/simx
'''

import sys, getopt
import pyCore

from mfem import path
from mfem.common.arg_parser import ArgParser 
import mfem.par as mfem
from mfem.par import intArray
from mfem.par import doubleArray
from mfem._par.pumi import ParPumiMesh
from mfem._par.pumi import ParMesh2ParPumiMesh
from os.path import expanduser, join
from mpi4py import MPI
import numpy as np
from numpy import sin, array, sqrt

model_file = ''
mesh_file  = ''
order = 1
geom_order = 1
adapt_ratio = 0.05
static_cond = False

freq = 1.0
kappa = np.pi * freq


class E_exact(mfem.VectorPyCoefficient):
    def __init__(self):
       mfem.VectorPyCoefficient.__init__(self, dim)
    def EvalValue(self, x):
       return (sin(kappa * x[1]),
               sin(kappa * x[2]),
               sin(kappa * x[0]))
class f_exact(mfem.VectorPyCoefficient):
    def __init__(self):
       mfem.VectorPyCoefficient.__init__(self, dim)
    def EvalValue(self, x):
       return ((1 + kappa**2)*sin(kappa * x[1]),
               (1 + kappa**2)*sin(kappa * x[2]),
               (1 + kappa**2)*sin(kappa * x[0]))

def get_field_z(pumi_mesh, field_name, field_type, template, grid):
  field_z = pyCore.createFieldOn(pumi_mesh, field_name, field_type)
  dim = pumi_mesh.getDimension()
  count_field = pyCore.createFieldOn(pumi_mesh, "count_field", pyCore.SCALAR)
  sol_field = pyCore.createFieldOn(pumi_mesh, "sol_field", pyCore.VECTOR)

  it = pumi_mesh.begin(0)
  while True:
    ent = pumi_mesh.iterate(it)
    if not ent:
      break
    pyCore.setScalar(count_field, ent, 0, 0.0)
    p = pyCore.Vector3(0.0, 0.0, 0.0)
    pyCore.setVector(sol_field, ent, 0, p)
  pumi_mesh.end(it)

  it = pumi_mesh.begin(dim)
  eid = 0
  while True:
    ent = pumi_mesh.iterate(it)
    if not ent:
      break

    vval = mfem.DenseMatrix()
    pmat = mfem.DenseMatrix()
    grid.GetVectorValues(eid, template, vval, pmat)
    pyCore.PCU_ALWAYS_ASSERT(vval.Width() == 4)
    pyCore.PCU_ALWAYS_ASSERT(vval.Height() == 3)

    for i in range(4):
      current_count = pumi_mesh.getVertScalarField(count_field, ent, i, 0)
      current_sol = pumi_mesh.getVertVectorField(sol_field, ent, i, 0)
      pumi_mesh.setVertScalarField(count_field, ent, i, 0, current_count + 1.0)
      pumi_mesh.setVertVectorField(sol_field, ent, i, 0, current_sol.x()+vval[0,i],
      	                                                 current_sol.y()+vval[1,i],
      	                                                 current_sol.z()+vval[2,i])
    eid = eid + 1

  pumi_mesh.end(it)

  it = pumi_mesh.begin(0)
  while True:
    ent = pumi_mesh.iterate(it)
    if not ent:
      break
    current_count = pyCore.getScalar(count_field, ent, 0)
    current_sol = pyCore.Vector3()
    pyCore.getVector(sol_field, ent, 0, current_sol)
    avg_sol = pyCore.Vector3(current_sol.x() / current_count,
    	                     current_sol.y() / current_count,
    	                     current_sol.z() / current_count)
    # mag = sqrt(avg_sol.x() * avg_sol.x() +
    # 	       avg_sol.y() * avg_sol.y() +
    # 	       avg_sol.z() * avg_sol.z())
    mag = avg_sol.x()
    pyCore.setScalar(field_z, ent, 0, mag)

  pumi_mesh.end(it)
  pumi_mesh.removeField(count_field)
  pyCore.destroyField(count_field)
  pumi_mesh.removeField(sol_field)
  pyCore.destroyField(sol_field)
  return field_z

def limit_refine_level(pumi_mesh, sizefield, level):
  it = pumi_mesh.begin(0)
  while True:
    ent = pumi_mesh.iterate(it)
    if not ent:
      break
    current_size = pumi_mesh.measureSize(ent)
    computed_size = pyCore.getScalar(sizefield, ent, 0)
    if computed_size < current_size / 2**level:
      computed_size = current_size / 2**level
    pyCore.setScalar(sizefield, ent, 0, computed_size)
  pumi_mesh.end(it)


parser = ArgParser(description='ex3p')
parser.add_argument('-m', '--mesh',
                    action = 'store', type = str,
                    help='Mesh file to use.')
parser.add_argument('-g', '--model',
                    action = 'store', type = str,
                    help='Model file to use.')
parser.add_argument('-o', '--order',
                    action = 'store', default = 1, type=int,
                    help = 'Finite Element Order')
parser.add_argument('-go', '--geom_order',
                    action = 'store', default = 1, type=int,
                    help = 'Geometric Representation Order')
parser.add_argument('-ar', '--adapt_ratio',
                    action = 'store', default = 0.05, type=float,
                    help = "adapt ratio")

args = parser.parse_args()
mesh_file = args.mesh
model_file = args.model
order = args.order
geom_order = args.geom_order
adapt_ratio = args.adapt_ratio

parser.print_options(args)

print 'Model     file is "', model_file
print 'Mesh      file is "', mesh_file
print 'order          is "', order
print 'geom_order     is "', geom_order
print 'adapt_ratio    is "', adapt_ratio


pyCore.lion_set_verbosity(1)
pyCore.PCU_Comm_Init()

num_proc = pyCore.PCU_Comm_Peers()
myid     = pyCore.PCU_Comm_Self()
verbose = (myid == 0)


pyCore.start_sim('simlog.txt')

pyCore.gmi_register_mesh()
pyCore.gmi_sim_start()
pyCore.gmi_register_sim()

pumi_mesh = pyCore.loadMdsMesh(model_file, mesh_file)
pyCore.printStats(pumi_mesh)
pumi_mesh.verify()

pyCore.writeASCIIVtkFiles('initial', pumi_mesh);

dim = pumi_mesh.getDimension()

pumi_mesh.verify()

pmesh = ParPumiMesh(pyCore.PCU_Get_Comm(), pumi_mesh) # supposed to me of type (ParMesh)
sdim = pmesh.SpaceDimension()

fec = mfem.ND_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)
size = fespace.GlobalTrueVSize()

if verbose: # note that size should be evaulated on all nodes
   print("Number of finite element unknowns: " + str(size))

ess_tdof_list = intArray();
if pmesh.bdr_attributes.Size():
    ess_bdr = intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);



b = mfem.ParLinearForm(fespace);
f = f_exact()
dd = mfem.VectorFEDomainLFIntegrator(f);
b.AddDomainIntegrator(dd)

x = mfem.ParGridFunction(fespace)
E = E_exact()
x.ProjectCoefficient(E);

muinv = mfem.ConstantCoefficient(1.0)
sigma = mfem.ConstantCoefficient(1.0)
a = mfem.ParBilinearForm(fespace)
a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(sigma))

if (static_cond):  a.EnableStaticCondensation()


ip_cnt = 0
vert_template = mfem.IntegrationRule(4)
ip = vert_template.IntPoint(0)
ip.Set3(0.0, 0.0, 0.0)
ip = vert_template.IntPoint(1)
ip.Set3(1.0, 0.0, 0.0)
ip = vert_template.IntPoint(2)
ip.Set3(0.0, 1.0, 0.0)
ip = vert_template.IntPoint(3)
ip.Set3(0.0, 0.0, 1.0)


i = 0
error = 10.0
while error > 0.25:
  if i > 10:
    break
  global_dofs = fespace.GlobalTrueVSize()
  if verbose:
    print "AMR Iteration: ", i
    print "Number of Unknowns: ", global_dofs

  a.Assemble()
  b.Assemble()

  ess_tdof_list = intArray();
  if pmesh.bdr_attributes.Size():
    ess_bdr = intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  A = mfem.HypreParMatrix()
  B = mfem.Vector()
  X = mfem.Vector()
  a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);


  prec_fespace = (a.SCParFESpace() if a.StaticCondensationIsEnabled() else fespace)
  ams = mfem.HypreAMS(A, prec_fespace)
  pcg = mfem.HyprePCG(A)
  pcg.SetTol(1e-12)
  pcg.SetMaxIter(500)
  pcg.SetPrintLevel(2)
  pcg.SetPreconditioner(ams)
  pcg.Mult(B, X)

  a.RecoverFEMSolution(X, b, x);

  err = x.ComputeL2Error(E)
  if verbose: # note that err should be evaulated on all nodes
      print("|| E_h - E ||_{L^2} = " + str(err))



  pmesh2 = ParMesh2ParPumiMesh(pmesh)

  field_z = get_field_z(pumi_mesh, "field_z", pyCore.SCALAR, vert_template, x)

  ipfield = pyCore.getGradIPField(field_z, "mfem_gradip", 2)
  sizefield = pyCore.getSPRSizeField(ipfield, adapt_ratio)

  limit_refine_level(pumi_mesh, sizefield, 1)

  pyCore.writeASCIIVtkFiles('before_adapt_' + str(i), pumi_mesh);

  pumi_mesh.removeField(ipfield)
  pyCore.destroyField(ipfield);

  pyCore.destroyNumbering(pumi_mesh.findNumbering("LocalVertexNumbering"));


  adapt_input = pyCore.configure(pumi_mesh, sizefield)
  adapt_input.shouldFixShape = True
  adapt_input.maximumIterations = 4
  pyCore.adapt(adapt_input)

  pyCore.writeASCIIVtkFiles('after_adapt_' + str(i), pumi_mesh);
  adapted_pmesh = ParPumiMesh(pyCore.PCU_Get_Comm(), pumi_mesh)
  pmesh2.UpdateMesh(adapted_pmesh)

  fespace.Update()
  x.Update()
  x.Assign(0)
  # pmesh2.FieldPUMItoMFEM(pumi_mesh, temp_field, x)

  x.ProjectCoefficient(E);
  a.Update()
  b.Update()

  pumi_mesh.removeField(field_z)
  pyCore.destroyField(field_z);
  pumi_mesh.removeField(sizefield)
  pyCore.destroyField(sizefield);
  i = i+1



smyid = '{:0>6d}'.format(myid)
mesh_name  =  "mesh."+smyid
sol_name   =  "ex6-sol."+smyid

# pmesh.PrintToFile(mesh_name, 8)
# x.SaveToFile(sol_name, 8)
