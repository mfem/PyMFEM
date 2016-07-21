%module ex3
%{
#include "linalg/vector.hpp"
%}

%callback("%s_cb");
%import "linalg/vector.hpp"
%inline %{
  
using namespace std;
using namespace mfem;

double freq = 1.0, kappa;
int dim;
 

void E_exact(const Vector &x, Vector &E)
{
   if (dim == 3)
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(2));
      E(2) = sin(kappa * x(0));
   }
   else
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(0));
      if (x.Size() == 3) { E(2) = 0.0; }
   }
}

void f_exact(const Vector &x, Vector &f)
{
   if (dim == 3)
   {
      f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
      f(1) = (1. + kappa * kappa) * sin(kappa * x(2));
      f(2) = (1. + kappa * kappa) * sin(kappa * x(0));
   }
   else
   {
      f(0) = (1. + kappa * kappa) * sin(kappa * x(1));
      f(1) = (1. + kappa * kappa) * sin(kappa * x(0));
      if (x.Size() == 3) { f(2) = 0.0; }
   }
}
%}
%nocallback;

%inline %{
void set_dim(const int in_dim)
{
   dim = in_dim;
}
 
void set_freq(const double in_freq)
{
   freq  = in_freq;
   kappa = freq * M_PI;
}
int get_dim(void)
{
   return dim;
}
double  get_freq(void)
{
   return freq;  
}
%}

