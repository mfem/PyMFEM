from __future__ import print_function

def nicePrint(*s):
    from mpi4py import MPI
    
    comm     = MPI.COMM_WORLD     
    nproc = MPI.COMM_WORLD.size
    myid     = MPI.COMM_WORLD.rank

    comm.Barrier() 
    for i in range(nproc):
        comm.Barrier()        
        if i == myid: print("[ID:"+str(myid)+"]" + ' '.join([str(ss) for ss in s]))
        comm.Barrier()
    
def niceCall(m, *args, **kargs):        
    from mpi4py import MPI
    
    comm     = MPI.COMM_WORLD     
    nproc = MPI.COMM_WORLD.size
    myid     = MPI.COMM_WORLD.rank

    comm.Barrier() 
    for i in range(nproc):
        comm.Barrier()        
        if i == myid:
            print("[ID:"+ str(myid) +"]")
            m(*args, **kargs)
        comm.Barrier()
        
