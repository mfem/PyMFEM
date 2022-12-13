/*
# 
#  using numba JIT function for mfem::FunctionCoefficient
#
(1) Using NumberFunction object  ---
(2) Using mfem.jit decorator (recommended)

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
#
#  This is a recommended way to use numba-compiled coefficient. It has
#     1) a simpler interface using mfem.jit decorator
#     2) dependecy support
#     3) complex function (a coefficient which returns complex)
#

# scalar coefficient
@mfem.jit.scalar(sdim=3)
def c12(ptx):
    return ptx[0] * ptx[sdim-1]  ### note sdim is defined when this is compiled

@mfem.jit.scalar(dependencies=((Er, Ei), density), complex=True))
def c12(ptx, E, density):
    return ptx[0] * (density * E.real + 1j*density.E.imag


# vectorr coefficient
@mfem.jit.vector(shape = (3,))
def f_exact(x, out):
    out_array = carray(out, (shape[0],))
    out[0] = (1 + kappa**2)*sin(kappa * x[1])
    out[1] = (1 + kappa**2)*sin(kappa * x[2])
    out[2] = (1 + kappa**2)*sin(kappa * x[0])

# passing Scalar/Vector/Matrix coefficient including GridFunctionCoefficients
# (Er, Ei) means complex number (GF for real and imaginary parts)
# density is double.

@mfem.jit.vector(dependencies=((Er, Ei), density), complex=True)
def f_exact(x, E, density, out):
    out[0] = (1 + kappa**2)*sin(kappa * x[1])
    out[1] = (1 + kappa**2)*sin(kappa * x[2])
    out[2] = (1 + kappa**2)*sin(kappa * x[0])

if return_complex=True, use .real and. imag as real and imaginary part
coefficient
   f_exact.real
   f_exact.imag
otherwise
   f_exact is coefficient

# vectorr coefficient
@mfem.jit.matrix(shape = (3,3))
def f_exact(x, out):
    out_array = carray(out, (shape[0], shape[1]))  
    for i in range(ndim):
        for j in range(ndim):
            out_array[i, j] = i*m + j

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
    double data[256];
    int datacount=0;
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
    double *GetData(){
      return data;
    }
    void SetDataCount(int x){
      return datacount=x;
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

//
//  FunctionCoefficeintExtra  (hold list of coefficients which is used as a parameter for function coefficient)
//
classe FunctionCoefficeintExtraBase:
 private:
  int num_coeffs;
  int num_vcoeffs;
  int num_mcoeffs;
  mfem::Coefficient *coeff;
  mfem::VectorCoefficient *vcoeff[];
  mfem::MatrizCoefficient *mcoeff[];

 protected:
  NumbaFunctionBase *obj);

 public:
  int ParamSize(){
    return size;
  }
  void SetParams(mfem::Coefficient *in_coeff[], int in_num_coeffs,
		 mfem::VectorCoefficient *in_vcoeff[], int in_num_vcoeffs,
		 mfem::MatrixCoefficient *in_mcoeff[], int in_num_mcoeffs){
    int size = 0;

    size += in_num_coeffs;
    for (int i = 0; i < in_num_vcoeffs; i++){
      size += vcoeff[i].GetVdim();
    }
    for (int i = 0; i < in_num_mcoeffs; i++){
      size += mcoeff[i].GetHeight() * mcoeff[i].GetWidth();
    }
    coeff = in_coeff;
    vcoeff = in_vcoeff;
    mcoeff = in_mcoeff;

    num_coeffs = in_num_coeffs;
    num_vcoeffs = in_num_vcoeffs;
    num_mcoeffs = in_num_mcoeffs;

    obj.SetDataCount(in_num_coeffs + in_num_vcoeffs + in_num_mcoeffs);

    MFEM_ASSERT(size < 256, "dependency dim must be less than 256");
  }

  double * PrepParams(){ElementTransformation &T,
                        const IntegrationPoint &ip){

    int vdim, h, w = 0;
    int idx = 0;

    double *data = obj.GetData();

    for (int i = 0; i < num_coeffs; i++){
      data[idx] = coeff[i].Eval(T, ip);
      idx++;
    }
    for (int i = 0; i < num_vcoeffs; i++){
      vdim = vcoeff[i].GetVdim();
      mfem::Vector V(vdim, &data[idx]);
      data[idx] = vcoeff[i].Eval(V, T, ip);
      idx += size;
    }
    for (int i = 0; i < num_mcoeffs; i++){
      h = mcoeff[i].GetHeight();
      w = mcoeff[i].GetWidth();
      mfem::DenseMatrix M(&data[idx], h, w);
      data[idx] = mcoeff[i].Eval(M, T, ip);
      idx += h*w;
    }
  }
};

classe FunctionCoefficeintExtra : public mfem::FunctionCoefficient,
  public FunctionCoefficeintExtra {
 public:
  FunctionCoefficeintExtra(std::function<double(const Vector &)> F,
			   NumbaFunctionBase *in_obj)
    : FunctionCoefficient(std::move(F)), obj(in_obj){}
 FunctionCoefficientExtra(std::function<double(const Vector &, double)> TDF,
			  NumbaFunctionBase *in_obj)
   : FunctionCoefficient(std::move(TDF)), obj(in_obj){}
  
 virtual double Eval(ElementTransformation &T,
		      const IntegrationPoint &ip){
   PrepParams(T, ip);
   return mfem::FunctionCoefficient::Eval(T, ip);
};

//  NumberFunction Implementation 2 (this is used for mfem.jit )
class NumbaFunction2 : public NumbaFunctionBase {
 private:
    std::function<double(const mfem::Vector &)> obj1;
    std::function<double(const mfem::Vector &, double t)> obj2;

 public:
    NumbaFunction(PyObject *input):
       NumbaFunctionBase(input, 3, false){}
    
    NumbaFunction(PyObject *input, bool td):
       NumbaFunctionBase(input, 3, td){}

    double call(const mfem::Vector &x){
      return ((double (*)(double *, double *))address_)(x.GetData(), data);
    }
    double callt(const mfem::Vector &x, double t){
      return ((double (*)(double *, double, double *))address_)(x.GetData(), t, data);
    }
    // complex return realpart
    double callr(const mfem::Vector &x){
      std::complex<double> ret;
      ret = ((std::complex<double> (*)(double *, double *))address_)(x.GetData(), data);
      return ret.real;
    }
    double calltr(const mfem::Vector &x, double t){
      std::complex<double> ret;            
      ret = ((std::complex<double> (*)(double *, double, double *))address_)(x.GetData(), t, data);
      return ret.real;      
    }
    // complex imag realpart
    double calli(const mfem::Vector &x){
      std::complex<double> ret;
      ret = ((std::complex<double> (*)(double *, double *))address_)(x.GetData(), data);
      return ret.imag;
    }
    double callti(const mfem::Vector &x, double t){
      std::complex<double> ret;            
      ret = ((std::complex<double> (*)(double *, double, double*, int))address_)(x.GetData(), t, data);
      return ret.imag;      
    }
    
    // FunctionCoefficient
    // mode   (0: real, 1: complex real part, 2: complex imag part)
    mfem::FunctionCoefficient* GenerateCoefficient(int mode=0){
      using std::placeholders::_1;
      using std::placeholders::_2;
      if (td_) {
	  switch(mode){
	  case 0:
	    obj2 = std::bind(&NumbaFunction::callt, this, _1, _2);
	  case 1:
	    obj2 = std::bind(&NumbaFunction::calltr, this, _1, _2);
	  case 2:
	    obj2 = std::bind(&NumbaFunction::callti, this, _1, _2);	    
          }
      } else {
	  switch(mode){
	  case 0:
	    obj2 = std::bind(&NumbaFunction::call, this, _1);
	  case 1:
	    obj2 = std::bind(&NumbaFunction::callr, this, _1);
	  case 2:
	    obj2 = std::bind(&NumbaFunction::calli, this, _1);	    
          }
      }    
      return new mfem::FunctionCoefficientExtra(obj1, this);
    }
};
 
 
%}

%newobject NumbaFunction::GenerateCoefficient;
%newobject VectorNumbaFunction::GenerateCoefficient;
%newobject MatrixNumbaFunction::GenerateCoefficient;
%newobject NumbaFunction2::GenerateCoefficient;
//%newobject VectorNumbaFunction::GenerateCoefficient;
//%newobject MatrixNumbaFunction::GenerateCoefficient;

%pythoncode %{

def generate_caller_scaler(settings):
    '''
    generate a callder function on the fly

    ex)
    if setting is {"input": (2, 1), "output": 2}

    def _caller(ptx, data):
        arr0 = data[0]+1j*data[0+1]
        arr2 = data[2]
        params = (arr0,arr2,)
        return (inner_func(ptx, *params))

    here inner_func is a function user provided.

    '''
    text = ['def _caller(ptx, data):']
    count = 0

    params_line = '    params = ('        
    for s in settings["input"]:
        if s == 2:
            t = '    arr'+str(count) + ' = data[' + str(count) + ']+1j*data[' + str(count) +'+1]'
            params_line += 'arr'+str(count)+','
            count = count + 2
        else:
            t = '    arr'+str(count) + ' = data[' + str(count) + ']'
            params_line += 'arr'+str(count)+','
            count = count + 1
        text.append(t)
    params_line += ')'

    text.append(params_line)
    text.append("    return (inner_func(ptx, *params))")
    return '\n'.join(text)
	      
def generate_signature_scalar(setting):
    '''
    generate a signature to numba-compile a user scalar function

    ex)
    if setting is {"input": (2, 1), "output": 2}
  
    output : types.complex128(CPointer(types.double), types.complex128,types.double,)

    '''

    sig = ''
    if setting['output'] == 1:
        sig += 'types.float64(CPointer(types.double, '
    else:
        sig += 'types.complex128(CPointer(types.double), '

    for s in setting['input']:
        if s == 1:
            sig += 'types.double,'
        else:
            sig += 'types.complex128,'

    sig = sig + ")"
    return sig

def generate_signature_array(setting):
    '''
    generate a signature to numba-compile a user vector/matrix function

    ex)
    if setting is {"input": (2, 1), "output": 2}
  
    output : types.void(CPointer(types.double), types.complex128, types.double, 
                        CPointer(types.complex128))

    '''

    sig = ''
    sig += 'types.void(CPointer(types.double, '
    for s in setting['input']:
        if s == 1:
            sig += 'types.double,'
        else:
            sig += 'types.complex128,'

    if setting['output'] == 1:
        sig += 'CPointer(types.double), '
    else:
        sig += 'CPointer(types.complex128), '


    sig = sig + ")"
    return sig
	      
def get_setting(iscomplex=False, dependencies=None):
    setting = {}
    if iscomplex:
        setting['output'] = 2
    else:
        setting['output'] = 1

    input = []
    for x in dependencis:
        if isinstance(x, tuple)
            input.append(2)
        else:
            input.append(1)

    setting["input"] = input
    return setting

try:
    from numba import cfunc, types, carray
    scalar_sig = types.double(types.CPointer(types.double),
			      types.intc,)
    scalar_sig_t = types.double(types.CPointer(types.double),
				types.double,
				types.intc,)
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

    class ComplexCoefficient()
        def __init__(self):
            self.real = None
	    self.imag = None

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

        def scalar(self, sdim=3, td=False, params={}, complex=False, dependencies=None):
            if dependencies is None:
                dependencies = []
            if params is None:
                params = {}
            params["sdim"] = sdim
			
            def dec(func):
	        from numba import cfunc, njit
			
                #l = len(signature(func).parameters) - len(dependencies)
	        setting = get_setting(complex, dependencies)

		sig = generate_signature_scalar(setting)
			
                gfunc=self._copy_func_and_apply_params(func, params)
                ff = ngit(sig)(gfunc)
			
                if td
                    caller_sig = types.double(types.CPointer(types.double),
				       types.double, 
				       types.CPointer(types.double))
		else:
                    caller_sig = types.double(types.CPointer(types.double),
                                        types.CPointer(types.double))

	        exec(generate_caller_scaler(settings))  # this defines _caller_func
	        callar_params = {'inner_func': ff}
                caller_func = self._copy_func_and_apply_params(_caller_func, caller_params)		      
                ff = cfunc(caller_sig)(caller_func)		      
                 

 	        if complex:
                     coeff = ComplexCoefficient()
	             coeff.real = NumbaFunction2(ff, td).GenerateCoefficient(1);
		     coeff.imag = NumbaFunction2(ff, td).GenerateCoefficient(2);		       
                else:
                     coeff = NumbaFunction(ff, sdim, td).GenerateCoefficient(0)
                return coeff
            return dec
       def vector(self, sdim=3, shape=None, td=False, params=None, complex=False, dependencies=None):
            shape = (sdim, ) if shape is None else shape
            if dependencies is None:
                dependencies = []
            if params is None:
                params = {}
            params["sdim"] = sdim
            params["shape"] = shape
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
      def matrix(self, sdim=3, shape=None, td=False, params=None, complex=False, dependencies=None):
            shape = (sdim, sdim) if shape is None else shape
            if dependencies is None:
                dependencies = []
            if params is None:
                params = {}
					   
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
