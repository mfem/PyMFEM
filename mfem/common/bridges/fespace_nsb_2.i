// Shared post-declaration bridge annotation for mfem::FiniteElementSpace.
NUMBA_SWIG_BRIDGE_CLASS(mfem::FiniteElementSpace)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/fespace/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/fespace/numba-swig-bridge.i"
#endif
#endif
