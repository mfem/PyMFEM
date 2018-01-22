import warnings
import numpy as np

from petram.mfem_config import use_parallel
if use_parallel:
    from petram.helper.mpi_recipes import *
    import mfem.par as mfem   
else:
    import mfem.ser as mfem

def make_matrix(x, y, z):
    pass

def do_findpoints(mesh, *args):
    sdim = mesh.SpaceDimension()
    
    shape = args[0].shape
    size= len(args[0].flatten())
    ptx = np.vstack([t.flatten() for t in args]).transpose().flatten()
    
    ptx2 = mfem.Vector(ptx)
    point_mat = mfem.DenseMatrix(ptx2.GetData(), size, sdim)
    elem_id = mfem.intArray()
    ips = mfem.IntegrationPointArray()
        
    num_found = mesh.FindPoints(point_mat, elem_id, ips, True)
    elem_id = np.array(elem_id.ToList())
    return v, elem_id, ips

def eval_at_points(gf, *args):
    args = [np.atleast_1d(np.array(t, copy=False)) for t in args]    
    mesh = gf.FESpace().Mesh()
    v, elem_id, ips = findpoints(mesh, *args)
    
def findpoints(mesh, *args):
    '''
    *args : x, y or x, y, z
    '''
    sdim = mesh.SpaceDimension()
    
    if len(args) != 3 and sdim == 3:
        assert False, "SpaceDimension = 3, pass x, y, z"
    elif len(args) != 2 and sdim == 2:
        assert False, "SpaceDimension = 3, pass x, y"
    elif len(args) != 1 and sdim == 1:
        assert False, "SpaceDimension = 3, pass x"
    else:
        args = [np.atleast_1d(np.array(t, copy=False)) for t in args]                
        return do_findpoints(mesh, *args)

