NUMBA_SWIG_BRIDGE_CLASS(mfem::CGSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::GMRESSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::IterativeSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::OperatorJacobiSmoother)
NUMBA_SWIG_BRIDGE_CLASS(mfem::OperatorChebyshevSmoother)
NUMBA_SWIG_BRIDGE_CLASS(mfem::SLISolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::FGMRESSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::BiCGSTABSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::MINRESSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::NewtonSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::LBFGSSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::OptimizationSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::SLBQPOptimizer)
NUMBA_SWIG_BRIDGE_CLASS(mfem::BlockILU)
NUMBA_SWIG_BRIDGE_CLASS(mfem::DirectSubBlockSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::ProductSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::OrthoSolver)
#ifdef MFEM_USE_SUITESPARSE
NUMBA_SWIG_BRIDGE_CLASS(mfem::UMFPackSolver)
NUMBA_SWIG_BRIDGE_CLASS(mfem::KLUSolver)
#endif
#ifdef MFEM_USE_MPI
NUMBA_SWIG_BRIDGE_CLASS(mfem::AuxSpaceSmoother)
#endif
#ifdef MFEM_USE_LAPACK
NUMBA_SWIG_BRIDGE_CLASS(mfem::NNLSSolver)
#endif
#ifdef NUMBA_SWIG_BRIDGE_GENERATED
#ifdef MFEM_NUMBA_SWIG_BRIDGE_PAR
%include "generated/mfem/_par/solvers/numba-swig-bridge.i"
#else
%include "generated/mfem/_ser/solvers/numba-swig-bridge.i"
#endif
#endif
