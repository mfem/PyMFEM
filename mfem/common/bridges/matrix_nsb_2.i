NUMBA_SWIG_BRIDGE_CLASS(mfem::Matrix)
NUMBA_SWIG_BRIDGE_CLASS(mfem::AbstractSparseMatrix)
#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/matrix/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/matrix/numba-swig-bridge.i"
#endif
#endif
