from mpi4py import MPI
import traceback
def get_mpi_datatype(data):
    '''
    get MPI datatype from numpy dtype
    '''
    #print '!!!chekcing MPI type ',  data.dtype.name
    s = data.dtype.itemsize  # get data size of element
    if data.dtype.name.startswith('int'):
        if s == 8:
           return MPI.LONG
        elif s == 4:
           return MPI.INT            
        elif s == 2:         
           return MPI.SHORT
        else:
           raise AttributeError("Unknow data type")
    if data.dtype.name.startswith('float'):       
        if s == 16:
           return MPI.LONG_DOUBLE
        elif s == 8:
           return MPI.DOUBLE
        else:
           raise AttributeError("Unknow data type")
       
