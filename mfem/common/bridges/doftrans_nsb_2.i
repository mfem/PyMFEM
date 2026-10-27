// Shared post-declaration bridge annotation for mfem::DofTransformation.
NUMBA_SWIG_BRIDGE_CLASS(mfem::DofTransformation)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/doftrans/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/doftrans/numba-swig-bridge.i"
#endif
#endif
