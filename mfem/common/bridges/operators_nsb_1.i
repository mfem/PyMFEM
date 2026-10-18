// Serial director annotation for the ex18 application operator.
%include "numba-swig-bridge-common.i"
%feature("nsb-director", "pymfem_TimeDependentOperator") mfem::TimeDependentOperator;
%feature("nsb-director", "pymfem_Operator") mfem::Operator;
%feature("nsb-director", "pymfem_Solver") mfem::Solver;
//NUMBA_SWIG_BRIDGE_DIRECTOR_OVERRIDE(
//    mfem::TimeDependentOperator::Mult)
