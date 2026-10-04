// Shared post-declaration bridge annotation for mfem::DenseMatrix.
NUMBA_SWIG_BRIDGE_CLASS(mfem::DenseMatrix)
NUMBA_SWIG_BRIDGE_CLASS(mfem::DenseTensor)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
%include "generated/mfem/_ser/densemat/numba-swig-bridge.i"
#endif
