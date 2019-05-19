#include "fem.hpp"
#include "coefficeint_py.hpp"

namespace mfem
{
  
void VectorFunctionCoefficientPy::Eval(Vector &V, ElementTransformation &T,
                                     const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   V.SetSize(vdim);
   if (Function)
   {
      EvalPy(transip, V);
   }
   else
   {
      EvalPyT(transip, GetTime(),  V);     

   }
   if (Q)
   {
      V *= Q->Eval(T, ip, GetTime());
   }
}
}  


