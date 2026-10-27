%include "numba-swig-bridge-common.i"
NUMBA_SWIG_BRIDGE_FUNCTION(mfem::Add)
NUMBA_SWIG_BRIDGE_FUNCTION(mfem::add)
NUMBA_SWIG_BRIDGE_FUNCTION(mfem::subtract)
NUMBA_SWIG_BRIDGE_FUNCTION(mfem::Mult)
NUMBA_SWIG_BRIDGE_IGNORE(mfem::Mult(const SparseMatrix&, const SparseMatrix&, SparseMatrix*))
