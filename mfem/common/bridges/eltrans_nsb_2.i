// Shared post-declaration bridge annotation for mfem::ElementTransformation.
NUMBA_SWIG_BRIDGE_CLASS(mfem::ElementTransformation)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/eltrans/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/eltrans/numba-swig-bridge.i"
#endif
#endif
