// Shared post-declaration bridge annotation for mfem::FiniteElement.
NUMBA_SWIG_BRIDGE_CLASS(mfem::FiniteElement)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/fe_base/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/fe_base/numba-swig-bridge.i"
#endif
#endif
