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


// NumbaCoefficient
/// NumbaFunction  : wrapper for numba compiled function
/// NumbaCoefficient : subclass of MFEM coefficient
class NumbaFunctionBase
{
 protected:
    void *address_;
    int  sdim_;
    bool td_;
    double data[256];
    double *data_ptr[16];
    int datacount=0;
 public:
  NumbaFunctionBase(PyObject *input, int sdim, bool td): sdim_(sdim), td_(td){ SetUserFunction(input); }
  void print_add(){ std::cout << address_ << '\n'; }
  void SetUserFunction(PyObject *input);
  double *GetData(){ return data; }
  double **GetPointer(){return data_ptr;}
  void SetDataCount(int x){datacount=x;}
  virtual ~NumbaFunctionBase(){}    
};

class NumbaCoefficientBase
{
 private:
  int num_coeffs;
  int num_vcoeffs;
  int num_mcoeffs;
  int num_dep = 0;
  int kinds[16];       // for now upto 16 coefficients
  int iscomplex[16];   // for now upto 16 coefficients
  mfem::Coefficient *coeff[16];
  mfem::VectorCoefficient *vcoeff[16];
  mfem::MatrixCoefficient *mcoeff[16];

 protected:
  NumbaFunctionBase *obj;

 public:
    NumbaCoefficientBase(NumbaFunctionBase *in_obj): obj(in_obj){}

  void SetParams(mfem::Coefficient *in_coeff[], int in_num_coeffs,
		 mfem::VectorCoefficient *in_vcoeff[], int in_num_vcoeffs,
		 mfem::MatrixCoefficient *in_mcoeff[], int in_num_mcoeffs);
  void PrepParams(mfem::ElementTransformation &T,
		  const mfem::IntegrationPoint &ip);
  void SetKinds(PyObject *kinds_);
  void SetIsComplex(PyObject *isComplex_);
  virtual ~NumbaCoefficientBase(){delete obj;}
};

class ScalarNumbaCoefficient : public mfem::FunctionCoefficient,  public NumbaCoefficientBase
{
 public:
  ScalarNumbaCoefficient(std::function<double(const mfem::Vector &)> F,
			   NumbaFunctionBase *in_obj)
    : FunctionCoefficient(std::move(F)), NumbaCoefficientBase(in_obj){}
  ScalarNumbaCoefficient(std::function<double(const mfem::Vector &, double)> TDF,
 			   NumbaFunctionBase *in_obj)
   : FunctionCoefficient(std::move(TDF)), NumbaCoefficientBase(in_obj){}
  
  virtual double Eval(mfem::ElementTransformation &T,
		      const mfem::IntegrationPoint &ip);
};

class VectorNumbaCoefficient : public mfem::VectorFunctionCoefficient,  public NumbaCoefficientBase
{
 public:
  VectorNumbaCoefficient(int dim, std::function<void(const mfem::Vector &, mfem::Vector &)> F,
		        NumbaFunctionBase *in_obj)
   : VectorFunctionCoefficient(dim, std::move(F)), NumbaCoefficientBase(in_obj){}
  VectorNumbaCoefficient(int dim, std::function<void(const mfem::Vector &, double, mfem::Vector &)> TDF,
   		        NumbaFunctionBase *in_obj)
    : VectorFunctionCoefficient(dim, std::move(TDF)), NumbaCoefficientBase(in_obj){}
  
 virtual void Eval(mfem::Vector &V,
		   mfem::ElementTransformation &T,
		   const mfem::IntegrationPoint &ip);
};

class MatrixNumbaCoefficient : public mfem::MatrixFunctionCoefficient,  public NumbaCoefficientBase
{
 public:
  MatrixNumbaCoefficient(int dim, std::function<void(const mfem::Vector &, mfem::DenseMatrix &)> F,
   			NumbaFunctionBase *in_obj)
   : MatrixFunctionCoefficient(dim, std::move(F)), NumbaCoefficientBase(in_obj){}
  MatrixNumbaCoefficient(int dim, std::function<void(const mfem::Vector &, double, mfem::DenseMatrix &)> TDF,
   			NumbaFunctionBase *in_obj)
    : MatrixFunctionCoefficient(dim, std::move(TDF)), NumbaCoefficientBase(in_obj){}
  
  virtual void Eval(mfem::DenseMatrix &K,
		    mfem::ElementTransformation &T,
		    const mfem::IntegrationPoint &ip);
};


