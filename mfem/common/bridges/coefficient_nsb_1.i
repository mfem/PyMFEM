// Direct public MFEM coefficient directors; references are callback-scoped.
%include "numba-swig-bridge-common.i"
%feature("nsb-director", "pymfem_Coefficient") mfem::Coefficient;
%feature("nsb-director", "pymfem_VectorCoefficient") mfem::VectorCoefficient;
%feature("nsb-director", "pymfem_MatrixCoefficient") mfem::MatrixCoefficient;
NUMBA_SWIG_BRIDGE_DIRECTOR_OVERRIDE(mfem::Coefficient::Eval(ElementTransformation &, const IntegrationPoint &))
NUMBA_SWIG_BRIDGE_DIRECTOR_OVERRIDE(mfem::VectorCoefficient::Eval(Vector &, ElementTransformation &, const IntegrationPoint &))
NUMBA_SWIG_BRIDGE_DIRECTOR_OVERRIDE(mfem::MatrixCoefficient::Eval(DenseMatrix &, ElementTransformation &, const IntegrationPoint &))
