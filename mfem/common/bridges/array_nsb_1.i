// Shared pre-declaration bridge setup for mfem::Array<int>.
%include "numba-swig-bridge-common.i"

NUMBA_SWIG_BRIDGE_ARRAY_VIEW_1D(
    mfem::Array<int>::GetDataArray,
    "dtype=int32;data=GetData;size=Size;mutable=true")
