#include <iostream>
#include "fem/intrules.hpp"

namespace mfem{
class PyIntegrationRule : public IntegrationRule
{
public:
  PyIntegrationRule():IntegrationRule(){}
  virtual ~PyIntegrationRule(){}
};
} // end of namespace
