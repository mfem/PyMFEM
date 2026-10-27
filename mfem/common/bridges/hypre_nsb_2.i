// Shared post-declaration bridge selections and generated adapters.
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreParVector)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreParMatrix)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreSmoother)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreTriSolve)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HyprePCG)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreGMRES)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreFGMRES)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreIdentity)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreDiagScale)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreParaSails)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreEuclid)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreILU)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreBoomerAMG)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreAMS)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreADS)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreLOBPCG)
NUMBA_SWIG_BRIDGE_CLASS(mfem::HypreAME)
#ifdef NUMBA_SWIG_BRIDGE_GENERATED
%include "generated/mfem/_par/hypre/numba-swig-bridge.i"
#endif
