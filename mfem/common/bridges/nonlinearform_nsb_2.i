// Shared post-declaration bridge annotation for mfem::NonlinearForm.
NUMBA_SWIG_BRIDGE_CLASS(mfem::NonlinearForm)
NUMBA_SWIG_BRIDGE_CLASS(mfem::BlockNonlinearForm)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/nonlinearform/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/nonlinearform/numba-swig-bridge.i"
#endif
#endif
