%module waveguide
%{
#include "linalg/vector.hpp"
%}

%callback("%s_cb");
%import "linalg/vector.hpp"
%inline %{
  
using namespace std;
using namespace mfem;

double freq = 1.0, kappa;
double zmid = 0.03;
int dim;
 

void Ht_port(const Vector &x, Vector &Ht)
{
     Ht(0) = 0.0;
     Ht(1) = 0.0;
     Ht(2) = cos(abs(x(2)-zmid)/zmid*M_PI_2);
}

void r_Jt_port(const Vector &x, Vector &Ht)
{
     Ht(0) = cos(abs(x(2)-zmid)/zmid*M_PI_2);
     Ht(1) = 0.0;
     Ht(2) = 0.0;
}

void i_Jt_port(const Vector &x, Vector &Ht)
{
     Ht(0) = 0.0;
     Ht(1) = 0.0;
     Ht(2) = 0.0;
}

void E_exact(const Vector &x, Vector &E)
{
   E(0) = 0.0;
   E(1) = 0.0;
   E(2) = 0.0;
}

%}
%nocallback;

%inline %{
void set_freq(const double in_freq)
{
   freq  = in_freq;
   kappa = freq * M_PI;
}
double  get_freq(void)
{
   return freq;  
}
%}

