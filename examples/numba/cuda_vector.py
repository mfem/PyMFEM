import mfem.ser as mfem
from numba import cuda
import numpy as np

@cuda.jit
def increment(values):
    i = cuda.grid(1)
    if i < values.size:
        values[i] += 1.0

data = np.arange(128*128, dtype=np.float64)

device = mfem.Device("cuda")
vector = mfem.Vector(data)

cuda.synchronize()  # Finish MFEM work before Numba accesses the storage.
values = vector.GetDeviceArray()
if values.size:
    increment[(values.size + 127) // 128, 128](values)
cuda.synchronize()  # Finish Numba work before returning to MFEM.

vector.HostRead()
print(vector.GetDataArray())  # [1. 2. 3. ... 128*128-1]