#
#   This example shows how to add waveguide port boundary 
#   condtion in plain MFEM and PyMFEM
#
#   This model simulates stragiht waveguide with 2 ports  
#   Bdry 2 : excitaiton port
#   Bdry 5 : passive port
#
#   Only TE01 modes are concerned, thought including multiple
#   higher modes is not complicated
#
#   In this module, Matrix is assembled in serial using serial
#   MFEM. Then, MUMPS is called in MPI mode
#  
#      DoF 0.47M (order = 2, refinement = 2), can be solved
#      with 4GB memory MacBookAir2014. Run time is about 1min.
#
#from ifigure.interactive import figure

#order 2  refine 2 = 472618
order = 2
ser_refine = 1

ext = 'port2'

from mfem import path    
import mfem.par as mfem
from mfem.par import intArray
from os.path import expanduser
import numpy as np
import scipy
import scipy.sparse
from mpi4py import MPI
from mfem.pi.utils import file_write, print_mem
from mfem.pi.solver.solver_utils import gather_vector, scatter_vector

epsilon0 = 8.8541878176e-12       # vacuum permittivity
mu0      = 4* np.pi*1e-7          # vacuum permiability
c        = 2.99792458e8           # speed of light
f        = 4.6e9                  # 5GHz
omega = 2*f*np.pi                 # omega =  2*pi*f
sigma = -0.0                  # resistivity
repsilon = 1.0                    # relative permittivity
rmu      = 1.0                    # relative permiability

m  = 0
n  = 1
a  = 0.007
b  = 0.06
k0 = 2*np.pi*f/c
k = np.sqrt(k0*k0 - np.pi*np.pi/b/b)
C = b/2
class C_Et(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi), 0, 0]
class C_Et2(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(2*x[2]/b*np.pi), 0, 0]
class C_Ht1(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [0., 0., k/mu0*np.sin(x[2]/b*np.pi)]
class C_Ht2(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [0., 0., -k/mu0*np.sin(x[2]/b*np.pi)]

num_proc = MPI.COMM_WORLD.size
myid     = MPI.COMM_WORLD.rank
if myid == 0:
    fid = open('waveguidep.out', "w")


meshfile = expanduser('~/src/PyMFEM/data/waveguide.mesh')
mesh = mfem.Mesh(meshfile, 1,1)
dim = mesh.Dimension()
sdim= mesh.SpaceDimension()

for x in range(ser_refine):
   mesh.UniformRefinement();
#mesh.ReorientTetMesh();

ess_tdof_list = intArray();
ess_bdr = intArray([1, 0, 1, 1, 0, 1])
                 
print_mem(myid)

pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
#del mesh                 
par_ref_levels = 0
for l in range(par_ref_levels):
    pmesh.UniformRefinement();
pmesh.ReorientTetMesh();
                 
fec = mfem.ND_FECollection(order, dim)
fespace = mfem.ParFiniteElementSpace(pmesh, fec)

fe_size = fespace.GlobalTrueVSize()
if (myid == 0):
   print('Number of finite element unknowns: '+  str(fe_size))

ess_tdof_list = mfem.intArray()
if pmesh.bdr_attributes.Size()>0:
#    ess_bdr = mfem.intArray(pmesh.bdr_attributes.Max())
#    ess_bdr.Assign(0)
    pass
fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

#from mfem.examples.waveguide import r_Jt_port_cb, i_Jt_port_cb, Ht_port_cb
r_b = mfem.ParLinearForm(fespace);
i_b = mfem.ParLinearForm(fespace);
port1_bdr = intArray([0, 1, 0, 0, 0, 0])
port2_bdr = intArray([0, 0, 0, 0, 1, 0])

#r_Jt = C_Jt(sdim)
if ext == 'port1':
   i_Ht = C_Ht1(sdim)
   i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port1_bdr)
else:
   i_Ht = C_Ht2(sdim)
   i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port2_bdr)
#r_dd = mfem.VectorFEBoundaryTangentLFIntegrator(r_f);
i_dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f); #\int test(E)*nxH dS
#r_b.AddBoundaryIntegrator(r_dd)
r_b.Assemble();
i_b.AddBoundaryIntegrator(i_dd)
i_b.Assemble();

    # 7. Define the solution vector x as a finite element grid function
#    corresponding to fespace. Initialize x by projecting the exact
#    solution. Note that only values from the boundary edges will be used
#    when eliminating the non-homogeneous boundary condition to modify the
#    r.h.s. vector b.

r_x = mfem.ParGridFunction(fespace)
i_x = mfem.ParGridFunction(fespace)

E = mfem.VectorArrayCoefficient(3)
for i in range(3):
    E.Set(i, mfem.ConstantCoefficient(0.0))
# Don't do this since, E will try to release the same object 3 times..    
#   zero = mfem.ConstantCoefficient(0.0)    
#   E.Set(0, zero);E.Set(1, zero);E.Set(2, zero)
class C_EexactR(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi)*np.cos(x[1]*k), 0, 0]
class C_EexactI(mfem.VectorPyCoefficient):
   def EvalValue(self, x):
       return [np.sin(x[2]/b*np.pi)*np.sin(x[1]*k), 0, 0]
ER = C_EexactR(sdim)
EI = C_EexactI(sdim)
# this takes a lot of time but needed only for verifing 
# the result with ideal solution. 
# r_x.ProjectCoefficient(ER);
# i_x.ProjectCoefficient(EI);
r_x.ProjectCoefficient(E);
i_x.ProjectCoefficient(E);
#(this is meaningless in parallel) sol_exact = r_x.GetDataArray() + 1j*i_x.GetDataArray()
                 
# 8. Set up the bilinear form corresponding to the EM diffusion operator
#       curl muinv curl + sigma I, by adding the curl-curl and the mass domain
#       integrators.
    
muinv = mfem.ConstantCoefficient(1./mu0/rmu);
    
mass = -(repsilon*epsilon0*omega*omega - 1j*omega*sigma)

co_rmass = mfem.ConstantCoefficient(mass.real)
co_imass = mfem.ConstantCoefficient(mass.imag)
    
r_a = mfem.ParBilinearForm(fespace);    # real part
i_a = mfem.ParBilinearForm(fespace);    # imag part
    
r_a.AddDomainIntegrator(mfem.CurlCurlIntegrator(muinv))
r_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(co_rmass))
i_a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(co_imass))
    
# 9. Assemble the bilinear form and the corresponding linear system,
#       applying any necessary transformations such as: eliminating boundary
#       conditions, applying conforming constraints for non-conforming AMR,
#       static condensation, etc.
    
r_a.Assemble()
i_a.Assemble()

print_mem(myid)
r_A = mfem.HypreParMatrix()
r_B = mfem.Vector()
r_X = mfem.Vector()
i_A = mfem.HypreParMatrix()
i_B = mfem.Vector()
i_X = mfem.Vector()

r_a.FormLinearSystem(ess_tdof_list, r_x, r_b, r_A, r_X, r_B)
i_a.FormLinearSystem(ess_tdof_list, i_x, i_b, i_A, i_X, i_B)

#r_A.Print('mat')
from mfem.pi.solver.mumps.hypre_to_mumps import sum_nnz, print_HYPRE_matrix_info
from mfem.pi.solver.mumps.hypre_to_mumps import ZMUMPS_LOC_Matrix, PyZMatrix, PyDMatrix
from mfem.pi.solver.mumps.hypre_to_mumps import get_true_local_nnz
from mfem.pi.solver.mumps.hypre_to_mumps import get_HypreParMatrixRow
from mfem.pi.solver.mumps.mumps_solve import z_real_array, d_array

BigM = ZMUMPS_LOC_Matrix(MPI.COMM_WORLD, 3,3)

#print get_true_local_nnz(r_A)
#print get_true_local_nnz(i_A)

BigM.add_real_hypre_matrix(r_A, 0)
BigM.add_imag_hypre_matrix(i_A, 0)

#for i in range(num_proc):
#    MPI.COMM_WORLD.Barrier()        
#    if ( myid == i ):  BigM.print_info()

#
# Port boundary
active_port = intArray([0, 1, 0, 0, 0, 0])
pt1 = mfem.ParLinearForm(fespace);
i_Ht = C_Ht1(sdim)
i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port1_bdr)
dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f)
pt1.AddBoundaryIntegrator(dd)
pt1.Assemble()
    
pt2 = mfem.ParLinearForm(fespace);
r_Et = C_Et(sdim)
r_f =  mfem.VectorRestrictedCoefficient(r_Et, port1_bdr)
dd = mfem.VectorFEDomainLFIntegrator(r_f)
pt2.AddBoundaryIntegrator(dd)
pt2.Assemble()
r_x1 = mfem.ParGridFunction(fespace)
r_x1.Assign(0.0)
r_x1.ProjectBdrCoefficientTangent(r_Et, port1_bdr)
#    pt2d = pt2.GetDataArray()/np.sum(pt2.GetDataArray()*r_x.GetDataArray())
    
passive_port = intArray([0, 0, 0, 0, 1, 0])
pt3 = mfem.ParLinearForm(fespace);
i_Ht = C_Ht2(sdim)
i_f =  mfem.VectorRestrictedCoefficient(i_Ht, port2_bdr)
dd = mfem.VectorFEBoundaryTangentLFIntegrator(i_f)
pt3.AddBoundaryIntegrator(dd)
pt3.Assemble()
    
pt4 = mfem.ParLinearForm(fespace);
r_Et = C_Et(sdim)
r_f =  mfem.VectorRestrictedCoefficient(r_Et, port2_bdr)
dd = mfem.VectorFEDomainLFIntegrator(r_f)
pt4.AddBoundaryIntegrator(dd)
pt4.Assemble()
r_x2 = mfem.ParGridFunction(fespace)
r_x2.Assign(0.0)
r_x2.ProjectBdrCoefficientTangent(r_Et, port2_bdr)
#    pt4d = pt4.GetDataArray()/np.sum(pt4.GetDataArray()*r_x.GetDataArray())

#for i in range(num_proc):
#    MPI.COMM_WORLD.Barrier()
#    if ( myid == i ): print myid, v.Size()
pt1b = pt1.ParallelAssemble().GetDataArray().copy()
pt2b = pt2.ParallelAssemble().GetDataArray().copy()
pt3b = pt3.ParallelAssemble().GetDataArray().copy()
pt4b = pt4.ParallelAssemble().GetDataArray().copy()
r_x1b= r_x1.ParallelAssemble().GetDataArray().copy()
r_x2b= r_x2.ParallelAssemble().GetDataArray().copy()
pt1d  = gather_vector(fespace, r_A, pt1b, MPI.DOUBLE, False)
pt2d  = gather_vector(fespace, r_A, pt2b, MPI.DOUBLE, False)
pt3d  = gather_vector(fespace, r_A, pt3b, MPI.DOUBLE, False)
pt4d  = gather_vector(fespace, r_A, pt4b, MPI.DOUBLE, False)
r_x1d  = gather_vector(fespace, r_A, r_x1b, MPI.DOUBLE, False)
r_x2d  = gather_vector(fespace, r_A, r_x2b, MPI.DOUBLE, False)
'''
for i in range(num_proc):
    MPI.COMM_WORLD.Barrier()
    if ( myid == i ): 
       print myid
       ess_tdof_list.Print()
       data = [get_HypreParMatrixRow(r_A, i) for i in ess_tdof_list.ToList()]
       print data
       data = [fespace.GetLocalTDofNumber(i) for i in ess_tdof_list.ToList()]
       print data
       data = [fespace.GetGlobalTDofNumber(i) for i in ess_tdof_list.ToList()]
       print data
       print len(data)
'''
data = np.array([get_HypreParMatrixRow(r_A, i) for i in ess_tdof_list.ToList()]).astype(np.int32)
all_ess = gather_vector(fespace, r_A, data, MPI.INT, False)
MPI.COMM_WORLD.Barrier()
if myid == 0:
    pt2d = pt2d/np.sum(pt2d*r_x1d)
    pt4d = pt4d/np.sum(pt4d*r_x2d)

    #print 'array check', np.sum(pt1d != 0.0), np.sum(pt2d != 0.0), np.sum(pt3d != 0.0), np.sum(pt4d != 0.0)

    pt1d= -pt1d
    pt3d = -pt3d
    #pt3d= pt3.GetDataArray()
    #print pt1d.shape, pt2d.shape, pt3d.shape, pt4d.shape
    #print list(all_ess)
    #print len(np.unique(all_ess))
    for x in list(all_ess):
       pt1d[x] = 0.0
       pt2d[x] = 0.0
       pt3d[x] = 0.0
       pt4d[x] = 0.0
    #print 'array check', np.sum(pt1d != 0.0), np.sum(pt2d != 0.0), np.sum(pt3d != 0.0), np.sum(pt4d != 0.0)
    pt1m = PyZMatrix(pt1d.shape[0], 1,  pt1d.shape[0])
    pt2m = PyZMatrix(1, pt2d.shape[0],  pt2d.shape[0])
    pt3m = PyZMatrix(pt3d.shape[0], 1,  pt3d.shape[0])
    pt4m= PyZMatrix(1, pt4d.shape[0],  pt4d.shape[0])
    
    pt1m.set_idata(z_real_array(pt1d))
    pt2m.set_rdata(z_real_array(pt2d))
    pt3m.set_idata(z_real_array(pt3d))
    pt4m.set_rdata(z_real_array(pt4d))

    #print pt1m.TrueNNZ(), pt2m.TrueNNZ(), pt3m.TrueNNZ(), pt4m.TrueNNZ()
    zero = PyZMatrix(1, 1, 0)
    one  = PyZMatrix(1, 1, 1)
    one.set_rdata(1.0)
    data = [pt1m, pt3m, pt2m, one, zero, pt4m, zero, one]
    for i, d in enumerate(data): BigM.add_py_matrix(d, i+1)

MPI.COMM_WORLD.Barrier()    
for i in range(8): BigM.share_py_matrix_info(i+1, 0)
BigM.assemble()

BigM.save_data('mat'+str(myid))
for i in range(num_proc):
    MPI.COMM_WORLD.Barrier()        
#    if ( myid == i ): BigM.print_info()
#if ( myid == 0 ): BigM.print_data()

# gather RHS

rhs_r = gather_vector(fespace, r_A, r_B.GetDataArray(), MPI.DOUBLE)
rhs_i = gather_vector(fespace, i_A, i_B.GetDataArray(), MPI.DOUBLE)

from mfem.pi.solver.mumps.mumps_solve import ZMUMPS, JOB_1_2_3, z_array, JOB_1_2
s = ZMUMPS()
s.set_icntl(5,0)
s.set_icntl(18,3)
s.set_icntl(28,2)
if myid ==0:
   s.set_n(BigM.N())
   rhs = rhs_r + 1j* rhs_i
   for x in list(all_ess):
       rhs[x] = 0.0
   if ext == 'port1':
      Bigb = np.hstack((rhs, 1.0, 0.0))
   else:
      Bigb = np.hstack((rhs, 0.0, 1.0))

   s.set_rhs(z_array(Bigb))
   #s.set_nz(BigM.NNZ())
   #s.set_irn(BigM.get_irn())
   #s.set_jcn(BigM.get_jcn())
   #s.set_a(BigM.get_data())

s.set_nz_loc(BigM.NNZ())
s.set_irn_loc(BigM.get_irn())
s.set_jcn_loc(BigM.get_jcn())
s.set_a_loc(BigM.get_data())

# No outputs
if False:
     s.set_icntl(1, -1)
     s.set_icntl(2, -1)
     s.set_icntl(3, -1)
     s.set_icntl(4,  0)

s.set_job(JOB_1_2_3)
s.run()
s.finish()

rsol = None; isol = None
if myid == 0:
   real_part = s.get_real_rhs()
   imag_part = s.get_imag_rhs()   
   rsol = real_part[:-2]
   isol = imag_part[:-2]
   ext = real_part[-2:] + 1j*imag_part[-2:]
   print 'R1, R2',  ext, np.abs(ext)
rdata = scatter_vector(fespace, rsol, MPI.DOUBLE)
idata = scatter_vector(fespace, isol, MPI.DOUBLE)

r_X.SetVector(mfem.Vector(rdata),0)
i_X.SetVector(mfem.Vector(idata),0)
r_a.RecoverFEMSolution(r_X, r_b, r_x)
i_a.RecoverFEMSolution(i_X, i_b, i_x)

smyid = '{:0>6d}'.format(myid)
mesh_name  =  "mesh."+smyid
rsol_name   =  "rsol."+smyid
isol_name   =  "isol."+smyid

pmesh.PrintToFile(mesh_name, 8)
r_x.SaveToFile(rsol_name, 8)
i_x.SaveToFile(isol_name, 8)
if 'ans' in locals(): ans(imag_part)

'''    
    from mfem.pi.solver.solver_utils import make_numpy_linear_system
    A1, b1 = make_numpy_linear_system(r_A, r_b)
    A2, b2 = make_numpy_linear_system(i_A, i_b)
    
    A = A1+1j*A2
    b = b1+1j*b2
    del A1, A2
    
    import scipy
    pt1d= pt1.GetDataArray()
    pt3d= pt3.GetDataArray()
    aaa = np.sum(pt2d*sol_exact)
    bbb = np.sum(pt4d*sol_exact)
    
    #v=figure();v.plot(pt1d,'r');v.plot(pt3d)
    sol_size = A.shape[0]
    BigA = scipy.sparse.bmat([[A, -1j*pt1d.reshape(-1,1), 1j*pt3d.reshape(-1,1)],
                              [pt2d.reshape(1,-1),  1,  0],
                              [pt4d.reshape(1,-1),  0,  1]], format='csr')
    del A
    BigExt = np.hstack((sol_exact, 0, bbb))
    Bigb = np.hstack((b, 1, 0))
    #v = figure();v.plot((BigA * BigExt).imag)
    #v = figure();v.plot((BigA * BigExt- Bigb).imag)
    ## Here, original version calls hegith, which is not
    ## defined in the header...!?
    file_write(fid, "Size of linear system: " , i_A.Size())

# 10. Solve
    BigA = BigA.tocoo(False)
    data = BigA.data.astype('complex')

    BigA.row+=1   #irn
    BigA.col+=1   #jcn
    rhs = Bigb
    #list(Bigb.astype('complex'))
    #if data.flags.contingous:
    #print BigA.row.flags
    #print BigA.col.flags
    #print BigA.data[:100]
    #print Bigb.flags
else:
   del mesh                 

MPI.COMM_WORLD.Barrier()
import mfem.pi.solver.mumps.mumps_solve as mumps_solve
s = mumps_solve.ZMUMPS()
s.set_icntl(28, 2)  ## parallel ordering
if myid ==0:

    file_write(fid, '!!!these two must be consistent')
    file_write(fid, 'sizeof(MUMPS_INT) ' , mumps_solve.SIZEOF_MUMPS_INT())
    file_write(fid, 'index data size ' , type(BigA.col[0]))
    file_write(fid, 'matrix data type ' , type(BigA.data[0]))
    file_write(fid, 'matrix size' , BigA.shape)
    s.set_n(BigA.shape[0])
    s.set_nz(len(data))
#    s.set_irn(mumps_solve.i_array(irn))
#    s.set_jcn(mumps_solve.i_array(jcn))
#    s.set_a(mumps_solve.z_array(list(data)))
    s.set_irn(mumps_solve.i_array(BigA.row))
    s.set_jcn(mumps_solve.i_array(BigA.col))
    s.set_a(mumps_solve.z_array(data))    
    s.set_rhs(mumps_solve.z_array(rhs))
     # No outputs
no_outputs = False
if no_outputs:
    s.set_icntl(1, -1)
    s.set_icntl(2, -1)
    s.set_icntl(3, -1)
    s.set_icntl(4,  0)
s.set_job(mumps_solve.JOB_1_2_3)
s.run()
s.finish()
if (myid == 0):
    #Bigsol = np.array(mumps_solve.z_to_list(s.get_rhs(), len(rhs)))
    Bigsol = Bigb
        
    sol = Bigsol[:sol_size]
    file_write(fid, 'R1, R2',  np.abs(BigExt[-2:]), np.abs(Bigsol[-2:]))
    r_X.SetVector(mfem.Vector(list(sol.real)), 0)
    i_X.SetVector(mfem.Vector(list(sol.imag)), 0)
    # 11. Recover the solution as a finite element grid function.
    r_a.RecoverFEMSolution(r_X, r_b, r_x)
    i_a.RecoverFEMSolution(i_X, i_b, i_x)
    
    mesh.PrintToFile('refined.mesh', 8)
    i_x.SaveToFile('soli.gf', 8)
    r_x.SaveToFile('solr.gf', 8)
    
    Exyz = [None]*3
    for idx in [1,2,3]:
       ReEx = mfem.Vector()
       r_x.GetNodalValues(ReEx, idx)  #index 1, 2, 3
       ImEx = mfem.Vector()
       i_x.GetNodalValues(ImEx, idx)  #index 1, 2, 3
       Exyz[idx-1] = ReEx.GetDataArray().copy() + 1j*ImEx.GetDataArray().copy()
'''
if (myid == 0):                 
    fid.close()

