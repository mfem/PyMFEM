#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
// common_functions uses a default-ignore allowlist; expose the resolver too.
%rename("%s") nsb_function_address_mfem__par_common_functions;
%include "generated/mfem/_par/common_functions/numba-swig-bridge.i"
#else
%rename("%s") nsb_function_address_mfem__ser_common_functions;
%include "generated/mfem/_ser/common_functions/numba-swig-bridge.i"
#endif
#endif
