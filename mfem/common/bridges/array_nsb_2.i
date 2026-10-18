// Shared post-declaration bridge annotation for mfem::Array<int>.
NUMBA_SWIG_BRIDGE_CLASS(mfem::Array<int>)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/array/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/array/numba-swig-bridge.i"
#endif
#endif
