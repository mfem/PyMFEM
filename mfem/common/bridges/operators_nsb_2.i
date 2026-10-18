// Include the generated serial director helpers after operator declarations.
NUMBA_SWIG_BRIDGE_CLASS(mfem::Operator)
#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/operators/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/operators/numba-swig-bridge.i"
#endif
#endif
