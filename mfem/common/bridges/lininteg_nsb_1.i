// Direct public MFEM linear-form integrator director.
%include "numba-swig-bridge-common.i"
%feature("nsb-director", "pymfem_LinearFormIntegrator") mfem::LinearFormIntegrator;
NUMBA_SWIG_BRIDGE_DIRECTOR_OVERRIDE(
    mfem::LinearFormIntegrator::AssembleRHSElementVect(
        const FiniteElement &, ElementTransformation &, Vector &))
