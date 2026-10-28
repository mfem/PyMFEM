// Python allocations are host memory regardless of MFEM's default allocator.
%pythonprepend mfem::Vector::Vector %{
from numpy import ndarray, ascontiguousarray
from .mem_manager import MemoryType
keep_link = False
own_data = False
if len(args) in (1, 2, 3) and isinstance(args[0], (list, ndarray)):
    if len(args) == 3 and args[2] != MemoryType.HOST:
        raise ValueError('Python storage requires MemoryType.HOST')
    source = args[0]
    if isinstance(source, ndarray):
        if source.dtype != 'float64':
            raise ValueError('Must be float64 array ' + str(source.dtype) + ' is given')
        if source.ndim != 1:
            raise ValueError('Vector input array must be one-dimensional')
        source = ascontiguousarray(source)
        keep_link = True
    else:
        own_data = True
    size = len(source) if len(args) == 1 else args[1]
    if size < 0 or size > len(source):
        raise ValueError('Vector size must fit the supplied host storage')
    args = (source, size, MemoryType.HOST)
%}

%pythonappend mfem::Vector::Vector %{
if keep_link:
    self._link_to_data = args[0]
if own_data:
    self.MakeDataOwner()
%}

%extend mfem::Vector {
    Vector(double *data_, int size_, mfem::MemoryType mt) {
        mfem::Memory<mfem::real_t> memory;
        memory.Wrap(data_, size_, mt, false);
        mfem::Vector *vector = new mfem::Vector();
        vector->NewMemoryAndSize(memory, size_, true);
        return vector;
    }
}
