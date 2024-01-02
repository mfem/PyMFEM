#
#   issue 3960/3965
#     Converting VisItDataCollection in serial format from an existing data collection saved in parallel
#     How do I Save a GridFunction created by GetSerialGridFunction
#
#   This code shows how to use GetSerialMehs and GetSerialGridFunction and to save the resultant
#   serial data in using VisItDataCollection
#
import os
import mfem.par as mfem
from mpi4py import MPI

num_procs = MPI.COMM_WORLD.size
myid = MPI.COMM_WORLD.rank
smyid = '{:0>6d}'.format(myid)

testDataFolder = os.path.abspath('./data_3960')

### import of pmesh
fname = mfem.MakeParFilename(testDataFolder+'/VisIt_volume/cycle_000000/pmesh.', myid)
pmesh = mfem.ParMesh(MPI.COMM_WORLD, fname)

### loading of the pressure field from the VisIt data collection
visit_dc = mfem.VisItDataCollection(testDataFolder+'/VisIt_volume/cycle')
visit_dc.SetMesh(MPI.COMM_WORLD, pmesh)
visit_dc.Load(0)
p = visit_dc.GetParField("pressure")

p.Save('ppressure.'+smyid)
pmesh.Print('pmesh.'+smyid)

smesh = pmesh.GetSerialMesh(0)
sp = p.GetSerialGridFunction(0, smesh);

if myid == 0:
   sp.Save('spressure')
   smesh.Print('smesh')
    
### printing additional information on the mesh
mesh_mfem = visit_dc.GetMesh()
mesh_data = mesh_mfem.GetNodes() # .GetTrueVector()  
order = int(mesh_mfem.GetNodalFESpace().FEColl().GetOrder())
dim = mesh_mfem.GetNodes().VectorDim()
dof = int(mesh_data.Size()/dim) # scalar degree of freedom, to be multiplied by number of scalar variables to have total dof
nele = mesh_mfem.GetNE() # .GetGlobalNE()
mesh_point_ordering = mesh_mfem.GetNodes().FESpace().GetOrdering() # "Ordering::byNODES=0" or "Ordering::byVDIM=1"

if myid == 0:
    print("dim="+str(dim))
    print("numb elem="+str(nele))
    print("order="+str(order))
    print("scalar dof="+str(int(dof)))
    print("pmesh_point_ordering ="+str(mesh_point_ordering))

### save the mesh in serial format
if myid == 0:
    folderList = [testDataFolder+'/VisIt_volume_serial/']
    for folderName in folderList:
        if not os.path.exists(folderName):
          os.mkdir(folderName)

if myid == 0:
    visit_dc_out = mfem.VisItDataCollection(testDataFolder+'/VisIt_volume_serial/cycle', smesh) 
    visit_dc_out.SetLevelsOfDetail(order)
    visit_dc_out.RegisterField("pressure", sp)
    visit_dc_out.SetCycle(0)
    visit_dc_out.SetTime(0)
    visit_dc_out.SetFormat(mfem.DataCollection.SERIAL_FORMAT)
    visit_dc_out.Save() 
