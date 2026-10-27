// Shared post-declaration bridge annotation for mfem::Vector.
NUMBA_SWIG_BRIDGE_CLASS(mfem::Vector)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/vector/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/vector/numba-swig-bridge.i"
#endif
#endif
