NUMBA_SWIG_BRIDGE_CLASS(mfem::Coefficient)
NUMBA_SWIG_BRIDGE_CLASS(mfem::VectorCoefficient)
NUMBA_SWIG_BRIDGE_CLASS(mfem::MatrixCoefficient)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/coefficient/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/coefficient/numba-swig-bridge.i"
#endif
#endif
