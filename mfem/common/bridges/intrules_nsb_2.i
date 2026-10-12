NUMBA_SWIG_BRIDGE_CLASS(mfem::IntegrationPoint)
NUMBA_SWIG_BRIDGE_CLASS(mfem::IntegrationRule)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/intrules/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/intrules/numba-swig-bridge.i"
#endif
#endif
