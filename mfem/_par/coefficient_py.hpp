#include "fem/coefficient.hpp"

class VectorFunctionCoefficientPy : public VectorFunctionCoefficient
{
private:


public:  
   using VectorFunctionCoefficient::Eval;
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);

   virtual void PyEval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);


}  
