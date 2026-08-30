'''
   MFEM example 2

   See c++ version in the MFEM library for more detail 
'''
import sys, getopt
import pyCore


# from mfem import path
import mfem.par as mfem
from mfem.par import intArray
from mfem.par import named_ifgzstream
from mfem.par import Vector
# from os.path import expanduser, join
# import numpy as np
from mfem._par.pumi import PumiMesh


def main(argv):
  model_file = ''
  mesh_file  = ''
  boundary_file  = ''
  try:
    opts, args = getopt.getopt(argv,"hg:m:b:",["model=","mesh=", "boundary="])
    print opts
    print args
  except getopt.GetoptError:
    print 'ex2.py -g <model> -m <mesh> -b <boundary>'
    sys.exit(3)
  for opt, arg in opts:
    print opt
    print arg
    if opt == '-h':
	print 'ex2.py -g <model> -m <mesh> -b <boundary'
	sys.exit()
    elif opt in ("-g", "--model"):
	model_file = arg
    elif opt in ("-m", "--mesh"):
	mesh_file = arg
    elif opt in ("-b", "--boundary"):
        boundary_file = arg
  print 'Model     file is "', model_file
  print 'Mesh      file is "', mesh_file
  print 'Boundary  file is "', boundary_file

  # LIONPRINT verbosity level to 1 for debugging
  pyCore.lion_set_verbosity(1)

  # PCU initialization
  pyCore.PCU_Comm_Init()

  # SIMX initialization
  pyCore.start_sim('simlog.txt')

  # gmi initialization
  pyCore.gmi_register_mesh()

  # gmi_sim start
  pyCore.gmi_sim_start()
  pyCore.gmi_register_sim()

  # load the pumi mesh and model and write the initial mesh to vtk
  pumi_mesh = pyCore.loadMdsMesh(model_file, mesh_file)
  pyCore.printStats(pumi_mesh)
  pumi_mesh.verify()
  pyCore.writeASCIIVtkFiles('before', pumi_mesh);

  # read boundary
  # for now set them manually.
  # TODO: figure out how to read the boundary file in python
  dirichlet = intArray()
  dirichlet.SetSize(2)
  dirichlet[0] = 1
  dirichlet[1] = 11
  dirichlet.Print()

  load = intArray()
  load.SetSize(1)
  load[0] = 5
  load.Print()



  # create a mfem mesh object form pumi mesh
  mfem_mesh = PumiMesh(pumi_mesh, 1, 1)
  dim = mfem_mesh.Dimension()

  print("HERE 02")
  it = pumi_mesh.begin(dim-1)
  bdr_cnt = 0
  while True:
    e = pumi_mesh.iterate(it)
    if not e: break
    model_type = pumi_mesh.getModelType(pumi_mesh.toModel(e))
    model_tag  = pumi_mesh.getModelTag(pumi_mesh.toModel(e))
    if model_type == (dim-1):
      mfem_mesh.GetBdrElement(bdr_cnt).SetAttribute(3)
      if dirichlet.Find(model_tag) != -1:
	mfem_mesh.GetBdrElement(bdr_cnt).SetAttribute(1)
      elif load.Find(model_tag) != -1:
	mfem_mesh.GetBdrElement(bdr_cnt).SetAttribute(2)
    bdr_cnt += 1
  pumi_mesh.end(it)

  for el in range(0, mfem_mesh.GetNE()):
    geom = mfem.Geometry()
    Tr = mfem_mesh.GetElementTransformation(el)
    ctr = Tr.Transform(geom.GetCenter(mfem_mesh.GetElementBaseGeometry(el)))
    if ctr[0] <= -0.05:
      mfem_mesh.SetAttribute(el, 1)
    elif ctr[0] >= 0.05:
      mfem_mesh.SetAttribute(el, 2)
    else:
      mfem_mesh.SetAttribute(el, 3)


  mfem_mesh.SetAttributes()
  if mfem_mesh.attributes.Max() < 2 or mfem_mesh.bdr_attributes.Max() < 2:
    sys.exit("Input mesh should have at leaset two materials and two boundary attributes!")

  print("HERE HERE")

  order = 1

  fec = mfem.H1_FECollection(order, dim)
  fespace = mfem.FiniteElementSpace(mfem_mesh, fec, dim)

  ess_tdof_list = intArray()
  ess_bdr = intArray([1]+[0]*(mfem_mesh.bdr_attributes.Max()-1))
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

  f = mfem.VectorArrayCoefficient(dim)
  for i in range(dim-1): f.Set(i, mfem.ConstantCoefficient(0.0));

  pull_force = mfem.Vector([0]*mfem_mesh.bdr_attributes.Max())
  pull_force[1] = -1.0e-2
  f.Set(dim-1, mfem.PWConstCoefficient(pull_force))
  f.Set(dim-2, mfem.PWConstCoefficient(pull_force))

  b = mfem.LinearForm(fespace)
  b.AddBoundaryIntegrator(mfem.VectorBoundaryLFIntegrator(f))
  b.Assemble()

  x = mfem.GridFunction(fespace)
  x.Assign(0.0)

  lamb = mfem.Vector(mfem_mesh.attributes.Max())
  lamb.Assign(1.0)
  lamb[0] = lamb[1]*10;
  lamb[1] = lamb[1]*100;
  lambda_func = mfem.PWConstCoefficient(lamb)

  mu = mfem.Vector(mfem_mesh.attributes.Max())
  mu.Assign(1.0);
  mu[0] = mu[1]*10;
  mu[1] = mu[1]*100;
  mu_func = mfem.PWConstCoefficient(mu)

  a = mfem.BilinearForm(fespace)
  a.AddDomainIntegrator(mfem.ElasticityIntegrator(lambda_func,mu_func))

  print('matrix...')
  static_cond = False
  if (static_cond): a.EnableStaticCondensation()
  a.Assemble()

  A = mfem.SparseMatrix()
  B = mfem.Vector()
  X = mfem.Vector()
  a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
  print('...done')## Here, original version calls hegith, which is not
  ## defined in the header...!?
  print("Size of linear system: " + str(A.Size()))

  M = mfem.GSSmoother(A)
  mfem.PCG(A, M, B, X, 1, 500, 1e-8, 0.0);

  a.RecoverFEMSolution(X, b, x);

  print("HERE --- ")

  if not mfem_mesh.NURBSext:
    mfem_mesh.SetNodalFESpace(fespace)

  nodes = mfem_mesh.GetNodes()
  nodes += x
  x *= -1
  mfem_mesh.PrintToFile('displaced.mesh', 8)
  x.SaveToFile('sol.gf', 8)


  # gmi_sim stop
  pyCore.gmi_sim_stop()

  # SIMX finalization
  pyCore.stop_sim()

  # gmi finalization
  pyCore.PCU_Comm_Free()


if __name__ == "__main__":
     main(sys.argv[1:])
