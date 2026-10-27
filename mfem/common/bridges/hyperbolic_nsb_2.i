// Shared post-declaration bridge annotations for ex18 hyperbolic objects.
NUMBA_SWIG_BRIDGE_CLASS(mfem::HyperbolicFormIntegrator)
NUMBA_SWIG_BRIDGE_CLASS(mfem::FluxFunction)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/hyperbolic/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/hyperbolic/numba-swig-bridge.i"
#endif
#endif
