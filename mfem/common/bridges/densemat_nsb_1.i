// Shared pre-declaration bridge setup for mfem::DenseMatrix.
%include "numba-swig-bridge-common.i"

NUMBA_SWIG_BRIDGE_ARRAY_VIEW_2D(
    mfem::DenseMatrix::GetDataArray,
    "dtype=float64;data=Data;rows=Height;cols=Width;layout=F;mutable=true")
NUMBA_SWIG_BRIDGE_ARRAY_VIEW_3D(
    mfem::DenseTensor::GetDataArray,
    "dtype=float64;data=Data;depth=SizeK;rows=SizeI;cols=SizeJ;layout=Kij;mutable=true")
