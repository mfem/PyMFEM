# CUDA arrays from MFEM Vectors

`vector.GetDeviceArray(access="readwrite")` returns a one-dimensional,
contiguous, float64 Numba CUDA `DeviceNDArray` sharing the Vector's device
storage. It also works on Vector-derived objects such as GridFunction and
ParticleVector. Serial and parallel bindings use the same API; each MPI rank
accesses its own local Vector.

This requires CUDA-enabled MFEM, an active CUDA backend, and Numba with working
CUDA support. Numba is imported only when the method is called. CPU, HIP,
single-precision bindings, and CUDA simulation are unsupported. MFEM and Numba
must use the same CUDA device/context. Standard `cuda` and `cuda:uvm` have been
validated; allocator integrations such as Umpire have not been tested here.

```python
import mfem.ser as mfem
from numba import cuda


@cuda.jit
def increment(values):
    i = cuda.grid(1)
    if i < values.size:
        values[i] += 1.0


device = mfem.Device("cuda")
vector = mfem.Vector([1.0, 2.0, 3.0])

cuda.synchronize()  # Finish MFEM work before Numba accesses the storage.
values = vector.GetDeviceArray()
if values.size:
    increment[(values.size + 127) // 128, 128](values)
cuda.synchronize()  # Finish Numba work before returning to MFEM.

vector.HostRead()
print(vector.GetDataArray())  # [2. 3. 4.]
```

The [particle update example](../examples/numba/device_array.py) advances all
local particles with one kernel launch per step:

```sh
python examples/numba/device_array.py
mpirun -n 2 python examples/numba/device_array.py -p
```

## Access and synchronization

| Access | MFEM acquisition | Caller responsibility |
| --- | --- | --- |
| `"read"` | `Read(True)` | Read existing values without modifying them. |
| `"readwrite"` (default) | `ReadWrite(True)` | Read or modify existing values. |
| `"write"` | `Write(True)` | Overwrite every element before reading or returning to MFEM. |

Read access advertises read-only metadata, but Numba currently does not enforce
it against kernel writes. The caller must honor the access mode.

The method does not synchronize CUDA work. Explicitly synchronize at both
handoffs, as above. The internal interface advertises the legacy default stream
and imports with `sync=False`; there is no custom-stream argument.

The returned view adds no separate data allocation or copy. Acquiring MFEM
device access can allocate device storage and transfer current host values.
Python list/NumPy constructors explicitly identify their storage as ordinary
host memory, including under `cuda:uvm`; such storage may require a separate
device allocation. MFEM-allocated Vectors retain the configured allocation
policy.

## Lifetime and memory validity

Acquire a **fresh view** after MFEM host access or any other change in memory
validity. Reusing a pointer does not notify MFEM of subsequent device writes.
Call `HostRead()` before reading results with `GetDataArray()`, or
`HostReadWrite()` before modifying them on the host.

The array retains the source Vector wrapper. For an owning Vector this keeps
the allocation alive through ordinary garbage collection. Borrowed Vectors and
subviews still require their external backing owners to remain alive.

Never resize, replace, free, or redistribute the backing storage while a view
is in use. Synchronize outstanding work, discard views, perform the operation,
and acquire new views afterward. The view does not lock the allocation.

An empty Vector produces a zero-length view; skip kernel launches for it.
Vector component ordering is preserved. For example, the particle example uses
component-major storage: all x values, followed by all y and all z values.

`GetDeviceArray` is called from Python before launching a CUDA kernel. It is
not a method callable inside Numba `@njit` functions or CUDA kernels, and it
requires no numba-swig-bridge registration.

## Validation commands

GPU tests are opt-in; the default invocation checks the CPU-side API only.

```sh
python test/test_device_array.py
python test/test_device_array.py --cuda
python test/test_device_array.py --managed
mpirun -n 2 python test/test_device_array.py -p --cuda
mpirun -n 2 python test/test_device_array.py -p --managed
python test/test_device_array.py --cuda --without-numba

# Constructor regressions independent of Numba:
python test/test_vector_host_memory.py --device cuda:uvm
mpirun -n 2 python test/test_vector_host_memory.py -p --device cuda:uvm
```
