#ifndef PYMFEM_PYBILININTEG
#define PYMFEM_PYBILININTEG

#include <iostream>
#include "fem/bilininteg.hpp"

namespace mfem{
class PyBilinearFormIntegrator : public BilinearFormIntegrator
{
public:
  PyBilinearFormIntegrator(const IntegrationRule *ir = NULL):BilinearFormIntegrator(ir){}
  virtual ~PyBilinearFormIntegrator() { }  
};
} // end of namespace

#endif //PYMFEM_PYBILININTEG
