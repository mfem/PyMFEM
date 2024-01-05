#include <iostream>
#include "fem/lininteg.hpp"

namespace mfem{
class PyLinearFormIntegrator : public LinearFormIntegrator
{
public:
  PyLinearFormIntegrator(const IntegrationRule *ir = NULL):LinearFormIntegrator(ir) {}    
  virtual ~PyLinearFormIntegrator() { }
};
} // end of namespace  

