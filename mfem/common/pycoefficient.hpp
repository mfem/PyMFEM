#include "fem/coefficient.hpp"

double fake_func(const mfem::Vector &);
void fake_func_vec(const mfem::Vector &, mfem::Vector &);
void fake_func_mat(const mfem::Vector &, mfem::DenseMatrix &);

namespace mfem{
class PyCoefficientBase : public FunctionCoefficient
{
 private:
   int isTimeDependent;
 public:
   /// Define a time-independent coefficient from a C-function
   PyCoefficientBase(int tdep): FunctionCoefficient(fake_func), isTimeDependent(tdep){}
   virtual void SetTime(double t){FunctionCoefficient::SetTime(t);}
   virtual double Eval(ElementTransformation &T,
                       const IntegrationPoint &ip);
   virtual double _EvalPy(Vector &)=0;
   virtual double _EvalPyT(Vector &, double)=0;
   virtual ~PyCoefficientBase() { }
};
class VectorPyCoefficientBase : public VectorFunctionCoefficient
{
private:
   int isTimeDependent;
public:
 VectorPyCoefficientBase(int dim,  int tdep, Coefficient *q=NULL): VectorFunctionCoefficient(dim, fake_func_vec, q), isTimeDependent(tdep) { }
   virtual void SetTime(double t){VectorFunctionCoefficient::SetTime(t);}  
   virtual void Eval(DenseMatrix &M, ElementTransformation &T,
		     const IntegrationRule &ir); /* do I need this method?? */
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);
   virtual void _EvalPy(Vector &, Vector &){};
   virtual void _EvalPyT(Vector &, double, Vector &){};
   virtual ~VectorPyCoefficientBase() { }
};
class MatrixPyCoefficientBase : public MatrixFunctionCoefficient{
private:
   int isTimeDependent;
public:
   MatrixPyCoefficientBase(int dim,  int tdep): MatrixFunctionCoefficient(dim, fake_func_mat), isTimeDependent(tdep) { }
   virtual void SetTime(double t){MatrixFunctionCoefficient::SetTime(t);}    
   virtual void Eval(DenseMatrix &K, ElementTransformation &T,
		     const IntegrationPoint &ip); 
   virtual void _EvalPy(Vector &, DenseMatrix &){};
   virtual void _EvalPyT(Vector &, double, DenseMatrix &){};
   virtual ~MatrixPyCoefficientBase() { }
};
}

