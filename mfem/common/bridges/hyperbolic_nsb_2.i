// Shared post-declaration bridge annotations for ex18 hyperbolic objects.
NUMBA_SWIG_BRIDGE_CLASS(mfem::HyperbolicFormIntegrator)
NUMBA_SWIG_BRIDGE_CLASS(mfem::FluxFunction)

#ifdef NUMBA_SWIG_BRIDGE_GENERATED
%include "generated/mfem/_ser/hyperbolic/numba-swig-bridge.i"
#endif
