'''
   MFEM example 6
      This is a version of Example 1 with a simple adaptive mesh
      refinement loop. 
      See c++ version in the MFEM library for more detail 
'''

import sys, getopt
import pyCore

from mfem import path
from mfem.common.arg_parser import ArgParser 
import mfem.par as mfem
from mfem.par import intArray
from mfem._par.pumi import ParPumiMesh
from mfem._par.pumi import ParMesh2ParPumiMesh
from os.path import expanduser, join
from mpi4py import MPI
import numpy as np

model_file = ''
mesh_file  = ''
order = 1
geom_order = 1
adapt_ratio = 0.05
static_cond = False


parser = ArgParser(description='ex6p')
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

pyCore.start_sim('simlog.txt')

pyCore.gmi_register_mesh()
pyCore.gmi_sim_start()
pyCore.gmi_register_sim()

pumi_mesh = pyCore.loadMdsMesh(model_file, mesh_file)
pyCore.printStats(pumi_mesh)
pumi_mesh.verify()

pyCore.writeASCIIVtkFiles('initial', pumi_mesh);

dim = pumi_mesh.getDimension()

ref_levels = int(np.floor(np.log(100000./pumi_mesh.count(dim))/np.log(2.)/dim))

if ref_levels > 2:
  uniform_refine = pyCore.configureUniformRefine(pumi_mesh, 1)
  pyCore.adapt(uniform_refine)

pumi_mesh.verify()

pyCore.writeASCIIVtkFiles('after_uniform_refine', pumi_mesh);

num_proc = pyCore.PCU_Comm_Peers()
myid     = pyCore.PCU_Comm_Self()
verbose = (myid == 0)

pmesh = ParPumiMesh(pyCore.PCU_Get_Comm(), pumi_mesh) # supposed to me of type (ParMesh)
sdim = pmesh.SpaceDimension()


fec = mfem.H1_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)

a = mfem.ParBilinearForm(fespace)
b = mfem.ParLinearForm(fespace)

one  = mfem.ConstantCoefficient(1.0)

integ = mfem.DiffusionIntegrator(one)
a.AddDomainIntegrator(integ)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(one))

x = mfem.ParGridFunction(fespace)
x.Assign(0)

if static_cond:
  a.EnableStaticCondensation()


max_iter = 3
for i in range(max_iter):
  global_dofs = fespace.GlobalTrueVSize()
  if myid == 0:
    print "AMR Iteration ", i
    print "Number of Unknowns ", global_dofs

  a.Assemble()
  b.Assemble()


  ess_tdof_list = intArray()
  if pmesh.bdr_attributes.Size() > 0:
    ess_bdr = intArray(pmesh.bdr_attributes.Max())
    ess_bdr.Assign(1)
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

  A = mfem.HypreParMatrix()
  B = mfem.Vector();  X = mfem.Vector()
  copy_interior = 1
  a.FormLinearSystem(ess_tdof_list, x, b, A, X, B, copy_interior)

  amg = mfem.HypreBoomerAMG()
  amg.SetPrintLevel(0)
  pcg = mfem.CGSolver(A.GetComm())
  pcg.SetPreconditioner(amg)
  pcg.SetOperator(A)
  pcg.SetRelTol(1e-6)
  pcg.SetMaxIter(200)
  pcg.SetPrintLevel(3)
  pcg.Mult(B, X)

  a.RecoverFEMSolution(X, b, x);

  tmag_field = pyCore.createFieldOn(pumi_mesh, "field_mag", pyCore.SCALAR)
  temp_field = pyCore.createFieldOn(pumi_mesh, "temp_field", pyCore.SCALAR)

  pmesh2 = ParMesh2ParPumiMesh(pmesh)
  pmesh2.FieldMFEMtoPUMI(pumi_mesh, x, temp_field, tmag_field)

  ipfield = pyCore.getGradIPField(tmag_field, "mfem_gradip", 2)
  sizefield = pyCore.getSPRSizeField(ipfield, adapt_ratio)


  pyCore.writeASCIIVtkFiles('before_adapt_' + str(i), pumi_mesh);

  pyCore.destroyField(ipfield);

  adapt_input = pyCore.configure(pumi_mesh, sizefield)
  adapt_input.shouldFixShape = True
  adapt_input.maximumIterations = 2
  pyCore.adapt(adapt_input)

  pyCore.writeASCIIVtkFiles('after_adapt_' + str(i), pumi_mesh);
  adapted_pmesh = ParPumiMesh(pyCore.PCU_Get_Comm(), pumi_mesh)
  pmesh2.UpdateMesh(adapted_pmesh)

  fespace.Update()
  x.Update()
  x.Assign(0)
  pmesh2.FieldPUMItoMFEM(pumi_mesh, temp_field, x)

  a.Update()
  b.Update()

  pyCore.destroyField(temp_field);
  pyCore.destroyField(tmag_field);
  pyCore.destroyField(sizefield);



smyid = '{:0>6d}'.format(myid)
mesh_name  =  "mesh."+smyid
sol_name   =  "ex6-sol."+smyid

pmesh.PrintToFile(mesh_name, 8)
x.SaveToFile(sol_name, 8)
