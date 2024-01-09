#ifndef PYMFEM_PYNONLININTEG
#define PYMFEM_PYNONLININTEG

#include <iostream>
#include "fem/nonlininteg.hpp"

namespace mfem{
class PyNonlinearFormIntegrator : public NonlinearFormIntegrator
{
public:
  PyNonlinearFormIntegrator(const IntegrationRule *ir = NULL):NonlinearFormIntegrator(ir) {}    
  virtual ~PyNonlinearFormIntegrator() { }
  
};
} // end of namespace

#endif //PYMFEM_PYNONLININTEG
