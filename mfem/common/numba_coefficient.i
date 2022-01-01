/*
# 
#  using numba JIT function for mfem::FunctionCoefficient
#
(1) Using NumberFunction object  ---
(2) Using mfem.jit decorator

--- (1) Using NumberFunction object  ---

from numba import cfunc, carray

from mfem.coefficient import (scalar_sig,
                              scalar_sig_t,
                              vector_sig,
                              vector_sig_t,
                              matrix_sig,
                              matrix_sig_t,

sdim = mesh.SpaceDimension()
vdim = mesh.Dimension()

@cfunc(scalar_sig)
def s_func(ptx, sdim):
    return 2*ptx[0]

c = NumbaFunction(s_func, sdim).GenerateCoefficient()

@cfunc(scalar_sig_t)
def s_func_t(ptx, t, sdim):
    return 2*ptx[0]
c = NumbaFunction(s_func, sdim, True).GenerateCoefficient()

@cfunc(vector_sig):
def v_func(ptx, out, sdim, vdim)
    out_array = carray(out, (vdim, ))
    for i in range(m):
        out_array[i] = x+i

c = VectorNumbaFunction(v_func, sdim, ndim).GenerateCoeficient()

@cfunc(vector_sig_t):
def v_func(ptx, t, out, sdim, vdim):
    out_array = carray(out, (vdim, ))
    for i in range(vdim):
        out_array[i] = x+i
c = VectorNumbaFunction(v_func, sdim, ndim, True).GenerateCoefficient()

@cfunc(matrix_sig):
def m_func(ptx, out, sdim, vdim):
    out_array = carray(out, (vdim, vdim))
    for i in range(ndim):
        for j in range(ndim):
            out_array[i, j] = i*m + j
c = MatrixNumbaFunction(m_func, sdim, ndim).GenerateCoefficient()

@cfunc(matrix_sig_t):
def m_func(ptx, t, out, sdim, vdim):
    out_array = carray(out, (vdim, vdim))
    for i in range(vdim):
        for j in range(vdim):
            out_array[i, j] = i*m + j
c = MatrixNumbaFunction(m_func, sdim, ndim, True).GenerateCoefficient()

--- (2) mfem.jit decorator ---
@mfem.jit.scalar()
def c12(ptx):
    return s_func0(ptx[0], ptx[1],  ptx[2])
@mfem.jit.vector()
def f_exact(x, out):
    out[0] = (1 + kappa**2)*sin(kappa * x[1])
    out[1] = (1 + kappa**2)*sin(kappa * x[2])
    out[2] = (1 + kappa**2)*sin(kappa * x[0])

Then, decorated function can be used as function coefficient
x.ProjectCoefficient(f_exact)       

*/

%pythonappend  NumbaFunction::GenerateCoefficient %{
    val._link = self
%}
%pythonappend  VectorNumbaFunction::GenerateCoefficient %{
    val._link = self
%}
%pythonappend  MatrixNumbaFunction::GenerateCoefficient %{
    val._link = self
%}

%inline %{
class NumbaFunctionBase{
 protected:
    void *address_;
    int  sdim_;
    bool td_;
 public:
    NumbaFunctionBase(PyObject *input, int sdim, bool td): sdim_(sdim), td_(td){
       PyObject* module = PyImport_ImportModule("numba.core.ccallback");
       if (!module){
           PyErr_SetString(PyExc_RuntimeError, "Can not load numba.core.ccallback module");
           return;
       }      
       PyObject* cls = PyObject_GetAttrString(module, "CFunc");
       if (!cls){
           PyErr_SetString(PyExc_RuntimeError, "Can not load CFunc");
           return;
       }
       int check = PyObject_IsInstance(input, cls);
       if (! check){
           PyErr_SetString(PyExc_RuntimeError, "Input must be numba.core.ccallback.CFunc");
           return;
       }
       PyObject *address = PyObject_GetAttrString(input, "address");
       void *ptr = PyLong_AsVoidPtr(address); 
       Py_DECREF(address);
       address_ = ptr;
    }
    double print_add(){
       std::cout << address_ << '\n';
       return 0;
    }
    virtual ~NumbaFunctionBase(){}    
};

class NumbaFunction : public NumbaFunctionBase {
 private:
    std::function<double(const mfem::Vector &)> obj1;
    std::function<double(const mfem::Vector &, double t)> obj2;

 public:
    NumbaFunction(PyObject *input, int sdim):
       NumbaFunctionBase(input, sdim, false){}
    
    NumbaFunction(PyObject *input, int sdim, bool td):
       NumbaFunctionBase(input, sdim, td){}

    double call0(const mfem::Vector &x){
      return ((double (*)(double *))address_)(x.GetData());
    }
    double call(const mfem::Vector &x){
      return ((double (*)(double *, int))address_)(x.GetData(), sdim_);
    }
    double call0t(const mfem::Vector &x, double t){
      return ((double (*)(double *, double))address_)(x.GetData(), t);
    }
    double callt(const mfem::Vector &x, double t){
      return ((double (*)(double *, double, int))address_)(x.GetData(), t, sdim_);
    }

    // FunctionCoefficient
    mfem::FunctionCoefficient* GenerateCoefficient(int use_0=0){
      using std::placeholders::_1;
      using std::placeholders::_2;
      if (use_0==0) {
	if (td_) {
            obj2 = std::bind(&NumbaFunction::callt, this, _1, _2);
            return new mfem::FunctionCoefficient(obj2);
        } else {
           obj1 = std::bind(&NumbaFunction::call, this, _1);
           return new mfem::FunctionCoefficient(obj1);
        }
      } else {
	if (td_) {
	  obj2 = std::bind(&NumbaFunction::call0t, this, _1, _2);
            return new mfem::FunctionCoefficient(obj2);
        } else {
	  obj1 = std::bind(&NumbaFunction::call0, this, _1);
           return new mfem::FunctionCoefficient(obj1);
        }
      }
   }
};
 
// VectorFunctionCoefficient     
class VectorNumbaFunction : public NumbaFunctionBase {
 private:
  std::function<void(const mfem::Vector &, mfem::Vector &)> obj1;
  std::function<void(const mfem::Vector &, double, mfem::Vector &)> obj2;
  int vdim_;
    
 public:
    VectorNumbaFunction(PyObject *input, int sdim, int vdim):
       NumbaFunctionBase(input, sdim, false), vdim_(vdim){}
    
    VectorNumbaFunction(PyObject *input, int sdim, int vdim, bool td):
       NumbaFunctionBase(input, sdim, td), vdim_(vdim){}

      
    void call(const mfem::Vector &x, mfem::Vector &out){
      out = 0.0;
      return ((void (*) (double *, double *, int, int))address_)(x.GetData(), 
    							      out.GetData(),
                                                              sdim_,
							      vdim_);      
       
    }
    void callt(const mfem::Vector &x, double t, mfem::Vector &out){
      out = 0.0;      
      return ((void (*) (double *, double,  double *, int, int))address_)(x.GetData(),
							              t,
								      out.GetData(),
								      sdim_,
	 							      vdim_);            
    }
    void call0(const mfem::Vector &x, mfem::Vector &out){
      out = 0.0;
      return ((void (*) (double *, double *))address_)(x.GetData(), 
						       out.GetData());
       
    }
    void call0t(const mfem::Vector &x, double t, mfem::Vector &out){
      out = 0.0;      
      return ((void (*) (double *, double,  double *))address_)(x.GetData(),
						                t,
								out.GetData());
    }

    mfem::VectorFunctionCoefficient* GenerateCoefficient(int use_0=0){
      using std::placeholders::_1;
      using std::placeholders::_2;
      using std::placeholders::_3;
      if (use_0 == 0){
         if (td_) {
    	   obj2 = std::bind(&VectorNumbaFunction::callt, this, _1, _2, _3);
	   return new mfem::VectorFunctionCoefficient(vdim_, obj2);
         } else {
    	   obj1 = std::bind(&VectorNumbaFunction::call, this, _1, _2);
	   return new mfem::VectorFunctionCoefficient(vdim_, obj1);
         }
      } else {
         if (td_) {	
    	   obj2 = std::bind(&VectorNumbaFunction::call0t, this, _1, _2, _3);
	   return new mfem::VectorFunctionCoefficient(vdim_, obj2);
         } else {
    	   obj1 = std::bind(&VectorNumbaFunction::call0, this, _1, _2);
	   return new mfem::VectorFunctionCoefficient(vdim_, obj1);
	 }
      }
   }
};
// MatrixFunctionCoefficient
class MatrixNumbaFunction : NumbaFunctionBase {
 private:
  std::function<void(const mfem::Vector &, mfem::DenseMatrix &)> obj1;
  std::function<void(const mfem::Vector &, double, mfem::DenseMatrix &)> obj2;
  int vdim_;
    
 public:
    MatrixNumbaFunction(PyObject *input, int sdim, int vdim):
       NumbaFunctionBase(input, sdim, false), vdim_(vdim){}
    MatrixNumbaFunction(PyObject *input, int sdim, int vdim, bool td):
       NumbaFunctionBase(input, sdim, td), vdim_(vdim){}

    void call(const mfem::Vector &x, mfem::DenseMatrix &out){
      out = 0.0;
      return ((void (*) (double *, double *, int, int))address_)(x.GetData(), 
    							      out.GetData(),
                                                              sdim_,
							      vdim_);      
       
    }
    void callt(const mfem::Vector &x, double t, mfem::DenseMatrix &out){
      out = 0.0;
      return ((void (*) (double *, double,  double *, int, int))address_)(x.GetData(),
							              t,
								      out.GetData(),
								      sdim_,
	 							      vdim_);            
    }
    void call0(const mfem::Vector &x, mfem::DenseMatrix &out){
      out = 0.0;
      return ((void (*) (double *, double *))address_)(x.GetData(), 
						       out.GetData());
       
    }
    void call0t(const mfem::Vector &x, double t, mfem::DenseMatrix &out){
      out = 0.0;
      return ((void (*) (double *, double,  double *))address_)(x.GetData(),
						                t,
							        out.GetData());
    }
       
    mfem::MatrixFunctionCoefficient* GenerateCoefficient(int use_0=0){
      using std::placeholders::_1;
      using std::placeholders::_2;
      using std::placeholders::_3;
      if (use_0 == 0) {
      if (td_) {
 	obj2 = std::bind(&MatrixNumbaFunction::callt, this, _1, _2, _3);
	return new mfem::MatrixFunctionCoefficient(vdim_, obj2);
      } else {
	obj1 = std::bind(&MatrixNumbaFunction::call, this, _1, _2);
	return new mfem::MatrixFunctionCoefficient(vdim_, obj1);
      }
      } else {
      if (td_) {
 	obj2 = std::bind(&MatrixNumbaFunction::call0t, this, _1, _2, _3);
	return new mfem::MatrixFunctionCoefficient(vdim_, obj2);
      } else {
	obj1 = std::bind(&MatrixNumbaFunction::call0, this, _1, _2);
	return new mfem::MatrixFunctionCoefficient(vdim_, obj1);
      }

      }
    }
};
 
%}

%newobject NumbaFunction::GenerateCoefficient;
%newobject NumbaFunction::GenerateCoefficientT;
%newobject VectorNumbaFunction::GenerateCoefficient;
%newobject VectorNumbaFunction::GenerateCoefficientT;
%newobject MatrixNumbaFunction::GenerateCoefficient;
%newobject MatrixNumbaFunction::GenerateCoefficientT;

%pythoncode %{
try:
    from numba import cfunc, types, carray
    scalar_sig = types.double(types.CPointer(types.double),
			  types.intc)
    scalar_sig_t = types.double(types.CPointer(types.double),
			    types.double,
			    types.intc)
    vector_sig = types.void(types.CPointer(types.double),
			types.CPointer(types.double),
			types.intc,
			types.intc)
    vector_sig_t = types.void(types.CPointer(types.double),
			  types.double,
                          types.CPointer(types.double),
			  types.intc,
			  types.intc)
  
    matrix_sig = vector_sig
    matrix_sig_t = vector_sig_t
      
    from inspect import signature
      
    class _JIT(object):
        def _copy_func_and_apply_params(self, f, params):
            import copy
            import types
            import functools
  
            """Based on https://stackoverflow.com/a/13503277/2988730 (@unutbu)"""
            globals = f.__globals__.copy()
            for k in params:
                globals[k] = params[k]
            g = types.FunctionType(f.__code__, globals, name=f.__name__,
                                   argdefs=f.__defaults__, closure=f.__closure__)
            g = functools.update_wrapper(g, f)
            g.__module__ = f.__module__
            g.__kwdefaults__ = copy.copy(f.__kwdefaults__)
            return g

        def func(self, sig, params):
            def dec(func):
                from numba import jit
                gfunc=self._copy_func_and_apply_params(func, params)
                ff = jit(sig)(gfunc)
                return ff
            return dec

        def scalar(self, sdim=3, td=False, params={}):
            def dec(func):
                l = len(signature(func).parameters)
                if l == 1 and not td:
                    sig = types.double(types.CPointer(types.double))
                    use_0 = 1
                elif l == 2 and not td:
                    sig = scalar_sig
                    use_0 = 0			  
                elif l == 2 and td:
                    sig = types.double(types.CPointer(types.double),
                    types.double)
                    use_0 = 1			  			  
                elif l == 3 and td:
                    sig = scalar_sig_t
                    use_0 = 0			  			  			  
                from numba import cfunc
                gfunc=self._copy_func_and_apply_params(func, params)
                ff = cfunc(sig)(gfunc)
                coeff = NumbaFunction(ff, sdim, td).GenerateCoefficient(use_0)
                return coeff
            return dec
        def vector(self, sdim=3, vdim=-1, td=False, params={}):
            vdim = sdim if vdim==-1 else vdim
            def dec(func):
                l = len(signature(func).parameters)
                if l == 2 and not td:
                    sig = types.void(types.CPointer(types.double),
                                     types.CPointer(types.double),)
                    use_0 = 1
                elif l == 4 and not td:
                    sig = vector_sig
                    use_0 = 0			  
                elif l == 3 and td:
                    sig = types.void(types.CPointer(types.double),
                                     types.CPointer(types.double),
                                     types.double)
                    use_0 = 1			  			  
                elif l == 5 and td:
                    sig = vector_sig_t
                    use_0 = 0			  			  			  
                else:
                    assert False, "Unsupported signature type"
                from numba import cfunc
                gfunc=self._copy_func_and_apply_params(func, params)		      
                ff = cfunc(sig)(gfunc)
                coeff = VectorNumbaFunction(ff, sdim, vdim, td).GenerateCoefficient(use_0)
                return coeff
            return dec
        def matrix(self, sdim=3, vdim=-1, td=False, params={}):
            vdim = sdim if vdim==-1 else vdim				   
            def dec(func):
                l = len(signature(func).parameters)
                if l == 2 and not td:
                    sig = types.void(types.CPointer(types.double),
                                     types.CPointer(types.double),)
                    use_0 = 1
                elif l == 4 and not td:
                    sig = matrix_sig
                    use_0 = 0			  
                elif l == 3 and td:
                    sig = types.void(types.CPointer(types.double),
                                     types.CPointer(types.double),
                                     types.double)
                    use_0 = 1			  			  
                elif l == 5 and td:
                    sig = matrix_sig_t
                    use_0 = 0			  			  			  
                else:
                    assert False, "Unsupported signature type"
                from numba import cfunc
                gfunc=self._copy_func_and_apply_params(func, params)
                ff = cfunc(sig)(gfunc)
                coeff = MatrixNumbaFunction(ff, sdim, vdim, td).GenerateCoefficient(use_0)
                return coeff
            return dec
    jit = _JIT()
except ImportError:
    pass
except BaseError:
    assert False, "Failed setting Numba signatures by an error other than ImportError"
%}
