// Shared pre-declaration bridge setup for mfem::Vector.
// This fragment is included immediately before linalg/vector.hpp.
%include "numba-swig-bridge-common.i"
%rename(_nsb_native_index) mfem::Vector::operator[];

NUMBA_SWIG_BRIDGE_ARRAY_VIEW_1D(
    mfem::Vector::GetDataArray,
    "dtype=float64;data=GetData;size=Size;mutable=true")
