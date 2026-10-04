// Serial director annotation for the ex18 application operator.
%include "numba-swig-bridge-common.i"
%feature("nsb-director", "pymfem_TimeDependentOperator") mfem::TimeDependentOperator;
NUMBA_SWIG_BRIDGE_DIRECTOR_OVERRIDE(
    mfem::TimeDependentOperator::Mult)
