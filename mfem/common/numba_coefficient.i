/*
#
#  USING numba JIT function for mfem::FunctionCoefficient
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
#     4) ptx and out are already processed using carray (no need to call it)
#
(usage)
@scalar(td=False, params={}, complex=False, dependency=None, interface='simple',
        sdim=None, debug=False)
@vector(vdim=None, shape=None, td=False, params={}, complex=False, dependency=None,
        interface='simple', sdim=None, debug=False)
@matrix(height=None, width=None, shape=None, td=False, params={}, complex=False, dependency=None,
        interface='simple', sdim=None, debug=False)

shape: shape of return value
vdim: vector dimensiton of vector coefficient
width and height: matrix size of matrix coefficient


td: time-dependence (False: stationary, True: time-dependent
complex: complex coefficient. if ture, it returns a tuple of coefficient
depenency: dependency to other coefficient
interface: calling proceture
   'simple': vector/matric function returns the result by value more pythonic.
   'c++': vector/matric function returns the result by parameter like C++
    other options: one can pass a tuple of (caller, signature) pair to create a custom
                   interface
sdim: space dimenstion (optional. if it is given, sdim can be used in the user function)
debug: extra debug print

(examples)
# scalar coefficient
@mfem.jit.scalar()
def c12(ptx):
    return ptx[0] * ptx[sdim-1]  ### note sdim is defined when this is compiled

@mfem.jit.scalar(dependency=((Er, Ei), density), complex=True))
def c12(ptx, E, density):
    return ptx[0] * (density * E[0].real + 1j*density.E[0].imag


# vectorr coefficient
@mfem.jit.vector(shape = (3,))
def f_exact(x, out):
    out[0] = (1 + kappa**2)*sin(kappa * x[1])
    out[1] = (1 + kappa**2)*sin(kappa * x[2])
    out[2] = (1 + kappa**2)*sin(kappa * x[0])

# passing Scalar/Vector/Matrix coefficient including GridFunctionCoefficients
# (Er, Ei) means complex number (GF for real and imaginary parts)
# density is double.

@mfem.jit.vector(dependency=((Er, Ei), density), complex=True)
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
@mfem.jit.matrix(3, shape = (3,3))
def f_exact(x, out):
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

%include "../common/typemap_macros.i"

LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::Coefficient*>&,
                         mfem::Coefficient*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::VectorCoefficient*>&,
                         mfem::VectorCoefficient*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<mfem::MatrixCoefficient*>&,
                         mfem::MatrixCoefficient*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<ScalarNumbaCoefficient*>&,
                         ScalarNumbaCoefficient*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<VectorNumbaCoefficient*>&,
                         VectorNumbaCoefficient*)
LIST_TO_MFEMOBJ_ARRAY_IN(const mfem::Array<MatrixNumbaCoefficient*>&,
                         MatrixNumbaCoefficient*)


%inline %{
void NumbaFunctionBase::SetUserFunction(PyObject *input){
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
//  NumbaCoefficientBase  (hold list of coefficients which is used as a parameter for function coefficient)
//
 template<>
   void NumbaCoefficientBase::SetParams(const mfem::Array<mfem::Coefficient *>& in_coeffs,
                                        const mfem::Array<mfem::VectorCoefficient *>& in_vcoeffs,
                                        const mfem::Array<mfem::MatrixCoefficient *>& in_mcoeffs,
                                        const mfem::Array<ScalarNumbaCoefficient *>&  in_nscoeffs,
                                        const mfem::Array<VectorNumbaCoefficient *>&  in_nvcoeffs,
                                        const mfem::Array<MatrixNumbaCoefficient *>&  in_nmcoeffs){

    static_assert(std::is_base_of<NumbaCoefficientBase, ScalarNumbaCoefficient>::value);
    static_assert(std::is_base_of<NumbaCoefficientBase, VectorNumbaCoefficient>::value);
    static_assert(std::is_base_of<NumbaCoefficientBase, MatrixNumbaCoefficient>::value);

    int size = 0;

    num_coeffs = in_coeffs.Size();
    num_vcoeffs = in_vcoeffs.Size();
    num_mcoeffs = in_mcoeffs.Size();

    if ((num_coeffs + num_vcoeffs + num_mcoeffs) > 16){
        throw std::invalid_argument("dependency dim must be up to 16");
    }
    if (num_mcoeffs > 16){
        throw std::invalid_argument("dependency dim must be up to 16");
    }
    if (num_vcoeffs > 16){
        throw std::invalid_argument("dependency dim must be up to 16");
    }
    if (num_coeffs > 16){
        throw std::invalid_argument("dependency dim must be up to 16");
    }

    pcoeffs = new mfem::Array<mfem::Coefficient *>(num_coeffs);
    pvcoeffs = new mfem::Array<mfem::VectorCoefficient *>(num_vcoeffs);
    pmcoeffs = new mfem::Array<mfem::MatrixCoefficient *>(num_mcoeffs);

    mfem::Array<mfem::Coefficient *>& coeffs = *pcoeffs;
    mfem::Array<mfem::VectorCoefficient *>& vcoeffs = *pvcoeffs;
    mfem::Array<mfem::MatrixCoefficient *>& mcoeffs = *pmcoeffs;

    for (int i = 0; i < num_coeffs; i++){
      coeffs[i] = in_coeffs[i];
      size ++;
    }
    for (int i = 0; i < num_vcoeffs; i++){
      vcoeffs[i] = in_vcoeffs[i];
      size += vcoeffs[i] -> GetVDim();
    }
    for (int i = 0; i < num_mcoeffs; i++){
      mcoeffs[i] = in_mcoeffs[i];
      size += mcoeffs[i] -> GetHeight() * mcoeffs[i] -> GetWidth();
    }

    obj->SetDataCount(num_coeffs + num_vcoeffs + num_mcoeffs);

    using std::invalid_argument;
    if (size > 256){
        throw std::invalid_argument("dependency dim must be less than 256");
    }

    num_ncoeffs = in_nscoeffs.Size();
    num_nvcoeffs = in_nvcoeffs.Size();
    num_nmcoeffs = in_nmcoeffs.Size();

    pncoeffs = new mfem::Array<NumbaCoefficientBase *>(num_ncoeffs);
    mfem::Array<NumbaCoefficientBase *>& ncoeffs = *pncoeffs;
    for (int i = 0; i < num_ncoeffs; i++){
      ncoeffs[i] = static_cast<NumbaCoefficientBase *>(in_nscoeffs[i]);
    }

    pnvcoeffs = new mfem::Array<NumbaCoefficientBase *>(num_nvcoeffs);
    mfem::Array<NumbaCoefficientBase *>& nvcoeffs = *pnvcoeffs;
    for (int i = 0; i < num_nvcoeffs; i++){
      nvcoeffs[i] = static_cast<NumbaCoefficientBase *>(in_nvcoeffs[i]);
    }

    pnmcoeffs = new mfem::Array<NumbaCoefficientBase *>(num_nmcoeffs);
    mfem::Array<NumbaCoefficientBase *>& nmcoeffs = *pnmcoeffs;
    for (int i = 0; i < num_nmcoeffs; i++){
      nmcoeffs[i] = static_cast<NumbaCoefficientBase *>(in_nmcoeffs[i]);
    }

 }

void NumbaCoefficientBase::SetTimeInDependency(double t){
    mfem::Array<mfem::Coefficient *>& coeffs = *pcoeffs;
    mfem::Array<mfem::VectorCoefficient *>& vcoeffs = *pvcoeffs;
    mfem::Array<mfem::MatrixCoefficient *>& mcoeffs = *pmcoeffs;

    for (int i = 0; i < num_coeffs; i++){
      coeffs[i]->SetTime(t);
    }
    for (int i = 0; i < num_vcoeffs; i++){
      vcoeffs[i]->SetTime(t);
    }
    for (int i = 0; i < num_mcoeffs; i++){
      mcoeffs[i]->SetTime(t);
    }

}
void NumbaCoefficientBase::PrepParams(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip){

    int vdim, h, w = 0;
    int idx = 0;
    double *data = obj -> GetData();
    double **data_ptr=obj ->GetPointer();

    int s_counter = 0;
    int v_counter = 0;
    int m_counter = 0;
    int ns_counter = 0;
    int nv_counter = 0;
    int nm_counter = 0;


    int counter = 0;

    mfem::Array<mfem::Coefficient *>& coeffs = *pcoeffs;
    mfem::Array<mfem::VectorCoefficient *>& vcoeffs = *pvcoeffs;
    mfem::Array<mfem::MatrixCoefficient *>& mcoeffs = *pmcoeffs;
    mfem::Array<NumbaCoefficientBase *>& ncoeffs = *pncoeffs;
    mfem::Array<NumbaCoefficientBase *>& nvcoeffs = *pnvcoeffs;
    mfem::Array<NumbaCoefficientBase *>& nmcoeffs = *pnmcoeffs;

    for (int i = 0; i < num_dep; i++){
        switch(kinds[i]){
        case 0:// scalar
          {
           data[idx] = coeffs[s_counter]->Eval(T, ip);
           data_ptr[counter] = &data[idx];

           idx ++;
           s_counter ++;
           counter ++;

           if (isdepcomplex[i] == 1){
               data[idx] = coeffs[s_counter]->Eval(T, ip);
               data_ptr[counter] = &data[idx];

               idx ++;
               s_counter ++;
               counter ++;
           }
           if (isdepcomplex[i] == 2){
             data[idx] = dynamic_cast<ScalarNumbaCoefficient *>(ncoeffs[ns_counter])->obj-> GetScalarImaginary();
              data_ptr[counter] = &data[idx];

              idx ++;
              ns_counter ++;
              counter ++;
           }
           break;
          }
        case 1:// vector
          {
           vdim = vcoeffs[i]->GetVDim();
           mfem::Vector V(vdim);
           vcoeffs[v_counter]->Eval(V, T, ip);

           data_ptr[counter] = &data[idx];
           for (int j = 0; j < vdim; j++){
             data[idx] =  V[j];
             idx ++;
           }
           v_counter ++;
           counter ++;

           if (isdepcomplex[i] == 1){
              vcoeffs[v_counter]->Eval(V, T, ip);

              data_ptr[counter] = &data[idx];
              for (int j = 0; j < vdim; j++){
                  data[idx] =  V[j];
                  idx ++;
              }
              v_counter ++;
              counter ++;
           }
           if (isdepcomplex[i] == 2){
              dynamic_cast<VectorNumbaCoefficient *>(nvcoeffs[nv_counter])->obj-> GetArrayImaginary(V);

              data_ptr[counter] = &data[idx];
              for (int j = 0; j < vdim; j++){
                  data[idx] =  V[j];
                  idx ++;
              }
              nv_counter ++;
              counter ++;
           }
           break;
          }
        case 2:// matrix
          {
           w = mcoeffs[m_counter]->GetWidth();
           h = mcoeffs[m_counter]->GetHeight();
           mfem::DenseMatrix M(h, w);
           mcoeffs[m_counter]->Eval(M, T, ip);

           data_ptr[counter] = &data[idx];
           for (int jj = 0; jj < w; jj++){
              for (int ii = 0; ii < h; ii++){
                 data[idx] =  M(ii, jj);
                 idx ++;
              }
           }
           m_counter ++;
           counter ++;

           if (isdepcomplex[i] == 1){
              mcoeffs[m_counter]->Eval(M, T, ip);
              data_ptr[counter] = &data[idx];
              for (int jj = 0; jj < w; jj++){
                for (int ii = 0; ii < h; ii++){
                  data[idx] =  M(ii, jj);
                  idx ++;
                }
              }
              m_counter ++;
              counter ++;
           }

           if (isdepcomplex[i] == 2){
              dynamic_cast<MatrixNumbaCoefficient *>(nmcoeffs[nm_counter])->obj-> GetArrayImaginary(M);

              data_ptr[counter] = &data[idx];
              for (int jj = 0; jj < w; jj++){
                for (int ii = 0; ii < h; ii++){
                  data[idx] =  M(ii, jj);
                  idx ++;
                }
              }
              nm_counter ++;
              counter ++;
           }

           break;
          }
        }
    }
}
void NumbaCoefficientBase::SetKinds(PyObject *kinds_){
  if (PyList_Check(kinds_)) {
     int ll = PyList_Size(kinds_);
     if (ll > 16){
       PyErr_SetString(PyExc_ValueError, "Dependecy must be less than 16");
       return;
     }
     for (int i = 0; i < ll; i++) {
        PyObject *s = PyList_GetItem(kinds_, i);
        kinds[i] = (int)PyInt_AsLong(s);
     }
     num_dep = ll;
  } else if (PyTuple_Check(kinds_)) {
     int ll = PyTuple_Size(kinds_);
     for (int i = 0; i < ll; i++) {
        PyObject *s = PyTuple_GetItem(kinds_,i);
        kinds[i] = (int)PyInt_AsLong(s);
     }
     if (ll > 16){
       PyErr_SetString(PyExc_ValueError, "Dependecy must be less than 16");
       return;
     }
     num_dep = ll;
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
  }
}
void NumbaCoefficientBase::SetIsDepComplex(PyObject *isComplex_){
  if (PyList_Check(isComplex_)) {
     int ll = PyList_Size(isComplex_);
     if (ll > 16){
       PyErr_SetString(PyExc_ValueError, "Dependecy must be less than 16");
       return;
     }
     for (int i = 0; i < ll; i++) {
        PyObject *s = PyList_GetItem(isComplex_, i);
        isdepcomplex[i] = (int)PyInt_AsLong(s);
     }
     num_dep = ll;
  } else if (PyTuple_Check(isComplex_)) {
     int ll = PyTuple_Size(isComplex_);
     for (int i = 0; i < ll; i++) {
        PyObject *s = PyTuple_GetItem(isComplex_,i);
        isdepcomplex[i] = (int)PyInt_AsLong(s);
     }
     num_dep = ll;
     if (ll > 16){
       PyErr_SetString(PyExc_ValueError, "Dependecy must be less than 16");
       return;
     }
  } else {
    PyErr_SetString(PyExc_ValueError, "Expecting a list/tuple");
  }
}


double ScalarNumbaCoefficient::Eval(mfem::ElementTransformation &T,
                                  const mfem::IntegrationPoint &ip){
   PrepParams(T, ip);
   return mfem::FunctionCoefficient::Eval(T, ip);
}

void VectorNumbaCoefficient::Eval(mfem::Vector &V,
                                  mfem::ElementTransformation &T,
                                  const mfem::IntegrationPoint &ip){
   V.SetSize(vdim);
   PrepParams(T, ip);
   return mfem::VectorFunctionCoefficient::Eval(V, T, ip);
  }

void MatrixNumbaCoefficient :: Eval(mfem::DenseMatrix &K,
                                    mfem::ElementTransformation &T,
                                    const mfem::IntegrationPoint &ip){
    K.SetSize(height, width);
    PrepParams(T, ip);
    return mfem::MatrixFunctionCoefficient::Eval(K, T, ip);
  }

//  NumberFunction Implementation 2 (this is used for mfem.jit )
class ScalarNumbaFunction2 : public NumbaFunctionBase {
 private:
  std::function<double(const mfem::Vector &)> obj1;
  std::function<double(const mfem::Vector &, double t)> obj2;
  std::complex<double> ret;

 public:
    ScalarNumbaFunction2(PyObject *input):
       NumbaFunctionBase(input, 3, false){}

    ScalarNumbaFunction2(PyObject *input, bool td):
       NumbaFunctionBase(input, 3, td){}

    ~ScalarNumbaFunction2(){}

    double call(const mfem::Vector &x){
      int sdim = x.Size();
      return ((double (*)(double *, int, void **))address_)(x.GetData(),
                                                            sdim,
                                                            (void **)data_ptr);
    }
    double callt(const mfem::Vector &x, double t){
      int sdim = x.Size();
      return ((double (*)(double *, int, double, void **))address_)(x.GetData(),
                                                                    sdim,
                                                                    t,
                                                                    (void **)data_ptr);
    }
    // complex real part
    double callr(const mfem::Vector &x){
      int sdim = x.Size();
      ret = ((std::complex<double> (*)(double *, int, void**))address_)(x.GetData(),
                                                                        sdim,
                                                                        (void **)data_ptr);
      return ret.real();
    }
    double calltr(const mfem::Vector &x, double t){
      int sdim = x.Size();
      ret = ((std::complex<double> (*)(double *, int, double, void**))address_)(x.GetData(),
                                                                                sdim,
                                                                                t,
                                                                                (void **)data_ptr);
      return ret.real();
    }
    // complex imag part
    double calli(const mfem::Vector &x){
      int sdim = x.Size();
      ret = ((std::complex<double> (*)(double *, int,  void**))address_)(x.GetData(),
                                                                         sdim,
                                                                         (void **)data_ptr);
      return ret.imag();
    }
    double callti(const mfem::Vector &x, double t){
      int sdim = x.Size();
      ret = ((std::complex<double> (*)(double *, int, double, void **))address_)(x.GetData(),
                                                                                 sdim,
                                                                                 t,
                                                                                 (void **)data_ptr);
      return ret.imag();
    }
    double GetScalarImaginary(){
      return ret.imag();
    }

    void set_obj1(std::function<double(const mfem::Vector &)> obj1_){
      obj1 = obj1_;
    };
    void set_obj2(std::function<double(const mfem::Vector &, double )> obj2_){
      obj2 = obj2_;
    };
    std::function<double(const mfem::Vector &)> get_obj1(){ return obj1; }
    std::function<double(const mfem::Vector &, double )> get_obj2(){return obj2; }
};

    // FunctionCoefficient
    // mode   (0: real, 1: complex real part, 2: complex imag part)
ScalarNumbaCoefficient* GenerateScalarNumbaCoefficient(PyObject *numba_func,  bool td, int mode){
      using std::placeholders::_1;
      using std::placeholders::_2;

      ScalarNumbaFunction2 *func_wrap = new ScalarNumbaFunction2(numba_func, td);
      if (td) {
          switch(mode){
          case 0:
            func_wrap->set_obj2(std::bind(&ScalarNumbaFunction2::callt, func_wrap, _1, _2));
            break;
          case 1:
            func_wrap->set_obj2(std::bind(&ScalarNumbaFunction2::calltr, func_wrap, _1, _2));
            break;
          case 2:
            func_wrap->set_obj2(std::bind(&ScalarNumbaFunction2::callti, func_wrap, _1, _2));
            break;
          }
          return new ScalarNumbaCoefficient(func_wrap->get_obj2(), func_wrap);
      } else {
          switch(mode){
          case 0:
            func_wrap->set_obj1(std::bind(&ScalarNumbaFunction2::call, func_wrap, _1));
            break;
          case 1:
            func_wrap->set_obj1(std::bind(&ScalarNumbaFunction2::callr, func_wrap, _1));
            break;
          case 2:
            func_wrap->set_obj1(std::bind(&ScalarNumbaFunction2::calli, func_wrap, _1));
            break;
          }
          return new ScalarNumbaCoefficient(func_wrap->get_obj1(), func_wrap);
      }
}
// VectorFunctionCoefficient
class VectorNumbaFunction2 : public NumbaFunctionBase {
 private:
  std::function<void(const mfem::Vector &, mfem::Vector &)> obj1;
  std::function<void(const mfem::Vector &, double, mfem::Vector &)> obj2;
  int vdim_;
  std::complex<double> *outc = nullptr;
 public:
    VectorNumbaFunction2(PyObject *input, int vdim)
       : NumbaFunctionBase(input, 3, false), vdim_(vdim){}

    VectorNumbaFunction2(PyObject *input, int vdim, bool td)
       : NumbaFunctionBase(input, 3, td), vdim_(vdim){}

    ~VectorNumbaFunction2(){
      delete [] outc;
    }

    void call(const mfem::Vector &x, mfem::Vector &out){
      int sdim = x.Size();
      out = 0.0;
      return ((void (*) (double *, int, void **, double *))address_)(x.GetData(),
                                                                     sdim,
                                                                     (void **)data_ptr,
                                                                     out.GetData());

    }
    void callt(const mfem::Vector &x, double t, mfem::Vector &out){
      int sdim = x.Size();
      out = 0.0;
      return ((void (*) (double *, int, double,  void**, double *))address_)(x.GetData(),
                                                                             sdim,
                                                                             t,
                                                                             (void **)data_ptr,
                                                                             out.GetData());
    }
    void callr(const mfem::Vector &x, mfem::Vector &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, void **, std::complex<double> *))address_)(x.GetData(),
                                                                            sdim,
                                                                            (void **)data_ptr, outc);
      for (int i = 0; i < vdim_; i++) {

        out[i] = outc[i].real();
      }
    }
    void calltr(const mfem::Vector &x, double t, mfem::Vector &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, double, void**,  std::complex<double> *))address_)(x.GetData(),
                                                                                    sdim,
                                                                                    t,
                                                                                    (void **)data_ptr,
                                                                                    outc);
      for (int i = 0; i < vdim_; i++) {
        out[i] = outc[i].real();
      }
    }
    void calli(const mfem::Vector &x, mfem::Vector &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, void**, std::complex<double> *))address_)(x.GetData(),
                                                                           sdim,
                                                                           (void **)data_ptr,
                                                                           outc);
      for (int i = 0; i < vdim_; i++) {
        out[i] = outc[i].imag();
      }

    }
    void callti(const mfem::Vector &x, double t, mfem::Vector &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, double, void**, std::complex<double> *))address_)(x.GetData(),
                                                                                   sdim,
                                                                                   t,
                                                                                   (void **)data_ptr,
                                                                                   outc);
      for (int i = 0; i < vdim_; i++) {
        out[i] = outc[i].imag();
      }
    }
    void GetArrayImaginary(mfem::Vector &out){
      for (int i = 0; i < vdim_; i++) {
        out[i] = outc[i].imag();
      }
    }

    void create_outc(){
      outc = new std::complex<double>[vdim_];
    }
    void set_obj1(std::function<void(const mfem::Vector &, mfem::Vector &)> obj1_){
      obj1 = obj1_;
    };
    void set_obj2(std::function<void(const mfem::Vector &, double, mfem::Vector &)> obj2_){
      obj2 = obj2_;
    };
    std::function<void(const mfem::Vector &, mfem::Vector &)> get_obj1(){return obj1;}
    std::function<void(const mfem::Vector &, double, mfem::Vector &)> get_obj2(){return obj2;}
};
VectorNumbaCoefficient* GenerateVectorNumbaCoefficient(PyObject *numba_func, int vdim, bool td, int mode){
      using std::placeholders::_1;
      using std::placeholders::_2;
      using std::placeholders::_3;

      VectorNumbaFunction2 *func_wrap = new VectorNumbaFunction2(numba_func, vdim, td);
      if (td) {
          switch(mode){
          case 0:
            func_wrap->set_obj2(std::bind(&VectorNumbaFunction2::callt, func_wrap, _1, _2, _3));
            break;
          case 1:
            func_wrap->set_obj2(std::bind(&VectorNumbaFunction2::calltr, func_wrap, _1, _2, _3));
            func_wrap->create_outc();
            break;
          case 2:
            func_wrap->set_obj2(std::bind(&VectorNumbaFunction2::callti, func_wrap, _1, _2, _3));
            func_wrap->create_outc();
            break;
          }
          return new VectorNumbaCoefficient(vdim, func_wrap->get_obj2(), func_wrap);
      } else {
          switch(mode){
          case 0:
            func_wrap->set_obj1(std::bind(&VectorNumbaFunction2::call, func_wrap, _1, _2));
            break;
          case 1:
            func_wrap->set_obj1(std::bind(&VectorNumbaFunction2::callr, func_wrap, _1, _2));
            func_wrap->create_outc();
            break;
          case 2:
            func_wrap->set_obj1(std::bind(&VectorNumbaFunction2::calli, func_wrap, _1, _2));
            func_wrap->create_outc();
            break;
          }
          return new VectorNumbaCoefficient(vdim, func_wrap->get_obj1(), func_wrap);
      }
}

// MatrixFunctionCoefficient
class MatrixNumbaFunction2 : public NumbaFunctionBase {
 private:
  std::function<void(const mfem::Vector &, mfem::DenseMatrix &)> obj1 = nullptr;
  std::function<void(const mfem::Vector &, double, mfem::DenseMatrix &)> obj2 = nullptr;
  int vdim_;
  std::complex<double> *outc;
 public:
    MatrixNumbaFunction2(PyObject *input, int vdim)
      : NumbaFunctionBase(input, 3, false), vdim_(vdim){}
    MatrixNumbaFunction2(PyObject *input, int vdim, bool td)
      : NumbaFunctionBase(input, 3, td), vdim_(vdim){}
    ~MatrixNumbaFunction2(){
      delete [] outc;
    }

    void call(const mfem::Vector &x, mfem::DenseMatrix &out){
      int sdim = x.Size();
      out = 0.0;
      return ((void (*) (double *, int, void**, double *))address_)(x.GetData(),
                                                                    sdim,
                                                                    (void **)data_ptr,
                                                                    out.GetData());

    }
    void callt(const mfem::Vector &x, double t, mfem::DenseMatrix &out){
      int sdim = x.Size();
      out = 0.0;
      return ((void (*) (double *, int, double, void**, double *))address_)(x.GetData(),
                                                                            sdim,
                                                                            t,
                                                                            (void **)data_ptr,
                                                                            out.GetData());
    }
    void callr(const mfem::Vector &x, mfem::DenseMatrix &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, void**, std::complex<double> *))address_)(x.GetData(),
                                                                           sdim,
                                                                           (void**)data_ptr,
                                                                           outc);
      double *outptr = out.GetData();
      for (int i = 0; i < vdim_; i++) {
        outptr[i] = outc[i].real();
      }
    }
    void calltr(const mfem::Vector &x, double t, mfem::DenseMatrix &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, double, void**,  std::complex<double> *))address_)(x.GetData(),
                                                                                    sdim,
                                                                                    t,
                                                                                    (void **)data_ptr,
                                                                                    outc);
      double *outptr = out.GetData();
      for (int i = 0; i < vdim_; i++) {
        outptr[i] = outc[i].real();
      }
    }
    void calli(const mfem::Vector &x, mfem::DenseMatrix &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, void**, std::complex<double> *))address_)(x.GetData(),
                                                                           sdim,
                                                                           (void **)data_ptr,
                                                                           outc);
      double *outptr = out.GetData();
      for (int i = 0; i < vdim_; i++) {
        outptr[i] = outc[i].imag();
      }

    }
    void callti(const mfem::Vector &x, double t, mfem::DenseMatrix &out){
      int sdim = x.Size();
      out = 0.0;
      ((void (*) (double *, int, double, void**,  std::complex<double> *))address_)(x.GetData(),
                                                                                    sdim,
                                                                                    t,
                                                                                    (void **)data_ptr,
                                                                                    outc);
      double *outptr = out.GetData();
      for (int i = 0; i < vdim_; i++) {
        outptr[i] = outc[i].imag();
      }
    }
    void GetArrayImaginary(mfem::DenseMatrix &out){
      double *outptr = out.GetData();
      for (int i = 0; i < vdim_; i++) {
        outptr[i] = outc[i].imag();
      }
    }
    void create_outc(){
      outc = new std::complex<double>[vdim_];
    }
    void set_obj1(std::function<void(const mfem::Vector &, mfem::DenseMatrix &)> obj1_){
      obj1 = obj1_;
    };
    void set_obj2(std::function<void(const mfem::Vector &, double, mfem::DenseMatrix &)> obj2_){
      obj2 = obj2_;
    };
    std::function<void(const mfem::Vector &, mfem::DenseMatrix &)> get_obj1(){return obj1;}
    std::function<void(const mfem::Vector &, double, mfem::DenseMatrix &)> get_obj2(){return obj2;}

};

MatrixNumbaCoefficient* GenerateMatrixNumbaCoefficient(PyObject *numba_func, int height,  int width, bool td, int mode){
  using std::placeholders::_1;
  using std::placeholders::_2;
  using std::placeholders::_3;

  //mfem::DenseMatrix &m = mfem::DenseMatrix(width, height);

  MatrixNumbaFunction2 *func_wrap = new MatrixNumbaFunction2(numba_func, width*height, td);
  if (td) {
          switch(mode){
          case 0:
            func_wrap->set_obj2(std::bind(&MatrixNumbaFunction2::callt, func_wrap, _1, _2, _3));
            break;
          case 1:
            func_wrap->set_obj2(std::bind(&MatrixNumbaFunction2::calltr, func_wrap, _1, _2, _3));
            func_wrap->create_outc();
            break;
          case 2:
            func_wrap->set_obj2(std::bind(&MatrixNumbaFunction2::callti, func_wrap, _1, _2, _3));
            func_wrap->create_outc();
            break;
          }
          return new MatrixNumbaCoefficient(height, width, func_wrap->get_obj2(), func_wrap);
   } else {
          switch(mode){
          case 0:
            func_wrap->set_obj1(std::bind(&MatrixNumbaFunction2::call, func_wrap, _1, _2));
            break;
          case 1:
            func_wrap->set_obj1(std::bind(&MatrixNumbaFunction2::callr, func_wrap, _1, _2));
            func_wrap->create_outc();
            break;
          case 2:
            func_wrap->set_obj1(std::bind(&MatrixNumbaFunction2::calli, func_wrap, _1, _2));
            func_wrap->create_outc();
            break;
          }
          return new MatrixNumbaCoefficient(height, width, func_wrap->get_obj1(), func_wrap);
      }
}

void SetNumbaCoefficientDependency(NumbaCoefficientBase &obj,
                                   const mfem::Array<mfem::Coefficient *>& a1,
                                   const mfem::Array<mfem::VectorCoefficient *>& a2,
                                   const mfem::Array<mfem::MatrixCoefficient *>& a3,
                                   const mfem::Array<ScalarNumbaCoefficient *>& a4,
                                   const mfem::Array<VectorNumbaCoefficient *>& a5,
                                   const mfem::Array<MatrixNumbaCoefficient *>& a6){
   obj.SetParams(a1,a2, a3, a4, a5, a6);
}
%}

%newobject NumbaFunction::GenerateCoefficient;
%newobject VectorNumbaFunction::GenerateCoefficient;
%newobject MatrixNumbaFunction::GenerateCoefficient;
%newobject GenerateScalarNumbaCoefficient;
%newobject GenerateVectorNumbaCoefficient;
%newobject GenerateMatrixNumbaCoefficient;

%pythoncode %{

from mfem.common.numba_coefficient_utils import (generate_caller_scalar,
                                                 generate_caller_array,
                                                 generate_caller_array_oldstyle,
                                                 generate_signature_scalar,
                                                 generate_signature_array,
                                                 generate_signature_array_oldstyle,
                                                 get_setting)


try:
    import numpy as np
    from numba import cfunc, types, carray, farray
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

    def IsNumbaCoefficient(obj):
        return (isinstance(obj, ScalarNumbaCoefficient) or
                isinstance(obj, VectorNumbaCoefficient) or
                isinstance(obj, MatrixNumbaCoefficient))

    def _copy_func_and_apply_params(f, params):
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

    class _JIT(object):
        def func(self, sig, params=None):
            params = {} if params is None else params
            def dec(func):
                from numba import jit
                gfunc=_copy_func_and_apply_params(func, params)
                ff = jit(sig)(gfunc)
                return ff
            return dec

        @staticmethod
        def scalar(func=None, **kwargs):
            def wrapper(func):
                def dec(*args, **kwargs):
                    return _scalar(*args, **kwargs)
                return dec(func, **kwargs)
            if func:
                return wrapper(func)
            else:
                return wrapper

        @staticmethod
        def vector(func=None, **kwargs):
            def wrapper(func):
                def dec(*args, **kwargs):
                    return _vector(*args, **kwargs)
                return dec(func, **kwargs)
            if func:
                return wrapper(func)
            else:
                return wrapper

        @staticmethod
        def matrix(func=None, **kwargs):
            def wrapper(func):
                def dec(*args, **kwargs):
                    return _matrixr(*args, **kwargs)
                return dec(func, **kwargs)
            if func:
                return wrapper(func)
            else:
                return wrapper

    jit = _JIT()
except ImportError:
    pass
except BaseError:
    assert False, "Failed setting Numba signatures by an error other than ImportError"

def _scalar(func, td=False, params=None, complex=False, dependency=None,
            interface="simple", sdim=None, debug=False):
    if dependency is None:
        dependency = []
    if params is None:
        params = {}
    params["sdim"] = sdim

    from numba import cfunc, njit

    setting = get_setting(1, complex, dependency, td)

    if interface=="c++":
        sig = generate_signature_scalar(setting)
    elif interface=="simple":
        sig = generate_signature_scalar(setting)
    else:
        sig = interface[1](setting)

    if debug:
        print("(DEBUG) signature for function:", sig)

    gfunc=_copy_func_and_apply_params(func, params)
    ff = njit(sig)(gfunc)

    if complex:
        outtype = types.complex128
    else:
        outtype = types.double

    if td:
        caller_sig = outtype(types.CPointer(types.double),
                             types.int32,
                             types.double,
                             types.CPointer(types.voidptr))
    else:
        caller_sig = outtype(types.CPointer(types.double),
                             types.int32,
                            types.CPointer(types.voidptr))

    if interface=="c++":
         caller_txt = generate_caller_scalar(setting)
    elif interface=="simple":
         caller_txt = generate_caller_scalar(setting)
    else:
      caller_txt = interface[0](setting)

    if debug:
         print("(DEBUG) generated caller function:\n", caller_txt)

    exec(caller_txt, globals(), locals())
    caller_params = {"inner_func": ff, "carray":carray, "farray":farray}
    caller_func = _copy_func_and_apply_params(locals()["_caller"], caller_params)
    ff = cfunc(caller_sig)(caller_func)

    if complex:
         coeff = GenerateScalarNumbaCoefficient(ff, td, 1)
         coeff.SetOutComplex(setting["output"])

         coeff.real = GenerateScalarNumbaCoefficient(ff, td, 1)
         coeff.imag = GenerateScalarNumbaCoefficient(ff, td, 2)
         coeffs = (coeff, coeff.real, coeff.imag)
    else:
         coeff = GenerateScalarNumbaCoefficient(ff, td, 0)
         coeff.SetOutComplex(setting["output"])
         coeffs = (coeff, )

    for c in coeffs:
         c.SetIsDepComplex(setting["isdepcomplex"])
         c.SetKinds(setting["kinds"])
         SetNumbaCoefficientDependency(c,
                                       setting["s_coeffs"],
                                       setting["v_coeffs"],
                                       setting["m_coeffs"],
                                       setting["ns_coeffs"],
                                       setting["nv_coeffs"],
                                       setting["nm_coeffs"])
         c._dependency_link = dependency
    return coeff


def _vector(func, vdim=None, shape=None, td=False, params=None,
            complex=False, dependency=None, interface="simple", sdim=None, debug=False):

    assert (vdim is not None or shape is not None), "vdim or shape must be given"

    if vdim is not None and shape is None:
        shape = (vdim, )
    if vdim is None and shape is not None:
        vdim = shape[0]

    assert vdim == shape[0], "vdim and shape are not consistent"

    if dependency is None:
       dependency = []
    if params is None:
       params = {}
    params["shape"] = shape
    params["vdim"] = vdim

    if sdim is not None:
       # optional parameter to use sdim in a user function
       params["sdim"] = sdim

    from numba import cfunc, njit

    setting = get_setting(shape, complex, dependency, td)

    if interface == "simple":
        sig = generate_signature_array(setting)
    elif interface == "c++":
        sig = generate_signature_array_oldstyle(setting)
    else:
        sig = interface[1](setting)

    if debug:
        print("(DEBUG) signature for function:", sig)

    gfunc=_copy_func_and_apply_params(func, params)
    ff = njit(sig)(gfunc)

    if complex:
        outtype = types.complex128
    else:
        outtype = types.double

    if td:
        caller_sig = types.void(types.CPointer(types.double),
                                types.int32,
                                types.double,
                                types.CPointer(types.voidptr),
                                types.CPointer(outtype))
    else:
        caller_sig = types.void(types.CPointer(types.double),
                                types.int32,
                                types.CPointer(types.voidptr),
                                types.CPointer(outtype))

    if interface == "simple":
        caller_txt = generate_caller_array(setting)
    elif interface == "c++":
        caller_txt = generate_caller_array_oldstyle(setting)
    else:
        caller_txt = interface[0](setting)

    if debug:
         print("(DEBUG) generated caller function:\n", caller_txt)

    exec(caller_txt, globals(), locals())

    caller_params = {"inner_func": ff, "np":np, "shape":shape,
                     "carray":carray, "farray":farray}

    if vdim is not None:
        caller_params["vdim"] = vdim

    caller_func = _copy_func_and_apply_params(locals()["_caller"], caller_params)
    ff = cfunc(caller_sig)(caller_func)

    if complex:
         coeff = GenerateVectorNumbaCoefficient(ff, shape[0], td, 1)
         coeff.SetOutComplex(setting["output"])

         coeff.real = GenerateVectorNumbaCoefficient(ff, shape[0], td, 1)
         coeff.imag = GenerateVectorNumbaCoefficient(ff, shape[0], td, 2)
         coeffs = (coeff, coeff.real, coeff.imag)
    else:
         coeff =  GenerateVectorNumbaCoefficient(ff, shape[0], td, 0)
         coeff.SetOutComplex(setting["output"])
         coeffs = (coeff, )

    for c in coeffs:
         c.SetIsDepComplex(setting["isdepcomplex"])
         c.SetKinds(setting["kinds"])
         SetNumbaCoefficientDependency(c,
                                       setting["s_coeffs"],
                                       setting["v_coeffs"],
                                       setting["m_coeffs"],
                                       setting["ns_coeffs"],
                                       setting["nv_coeffs"],
                                       setting["nm_coeffs"])
         c._dependency_link = dependency
    return coeff


def _matrix(func, height=None, width=None, shape=None, td=False, params=None,
            complex=False, dependency=None, interface="simple", sdim=None, debug=False):

    if (width is None and height is not None or
        width is not None and height is None) :
        assert False, "height and width must be used together"

    assert (width is not None or shape is not None), "w/h or shape must be given"

    if width is not None and shape is None:
        shape = (width, height)
    if height is None and width is None and shape is not None:
        width = shape[0]
        height = shape[1]
    assert height == shape[0], "height and shape[0] are not consistent"
    assert width == shape[1], "width and shape[1] are not consistent"

    if shape[0] != shape[1]:
        import warnings
        warning.warn("Rectangular matrix coefficient is experimental", UserWarning)

    if dependency is None:
        dependency = []
    if params is None:
        params = {}
    params["sdim"] = sdim
    params["shape"] = shape
    params["width"] = width
    params["height"] = height

    from numba import cfunc, njit

    setting = get_setting(shape, complex, dependency, td)

    if interface == "simple":
        sig = generate_signature_array(setting)
    elif interface == "c++":
        sig = generate_signature_array_oldstyle(setting)
    else:
        sig = interface[1](setting)

    if debug:
        print("(DEBUG) signature for function:", sig)

    gfunc = _copy_func_and_apply_params(func, params)
    ff = njit(sig)(gfunc)

    if complex:
        outtype = types.complex128
    else:
        outtype = types.double
    if td:
        caller_sig = types.void(types.CPointer(types.double),
                                types.int32,
                                types.double,
                                types.CPointer(types.voidptr),
                                types.CPointer(outtype))
    else:
        caller_sig = types.void(types.CPointer(types.double),
                                types.int32,
                                types.CPointer(types.voidptr),
                                types.CPointer(outtype))

    if interface == "simple":
        caller_txt = generate_caller_array(setting)
    elif interface == "c++":
        caller_txt = generate_caller_array_oldstyle(setting)
    else:
        caller_txt = interface[0](setting)

    if debug:
         print("(DEBUG) generated caller function:\n", caller_txt)

    exec(caller_txt, globals(), locals())

    caller_params = {"inner_func": ff, "np":np, "shape":shape,
                     "carray":carray, "farray":farray}
    caller_func = _copy_func_and_apply_params(locals()["_caller"], caller_params)
    ff = cfunc(caller_sig)(caller_func)

    if complex:
         coeff = GenerateMatrixNumbaCoefficient(ff, shape[0], shape[1], td, 1)
         coeff.SetOutComplex(setting["output"])

         coeff.real = GenerateMatrixNumbaCoefficient(ff, shape[0], shape[1], td, 1)
         coeff.imag = GenerateMatrixNumbaCoefficient(ff, shape[0], shape[1], td, 2)
         coeffs = (coeff, coeff.real, coeff.imag)

    else:
         coeff = GenerateMatrixNumbaCoefficient(ff, shape[0], shape[1], td, 0)
         coeff.SetOutComplex(setting["output"])
         coeffs = (coeff, )

    for c in coeffs:
         c.SetIsDepComplex(setting["isdepcomplex"])
         c.SetKinds(setting["kinds"])
         SetNumbaCoefficientDependency(c,
                                       setting["s_coeffs"],
                                       setting["v_coeffs"],
                                       setting["m_coeffs"],
                                       setting["ns_coeffs"],
                                       setting["nv_coeffs"],
                                       setting["nm_coeffs"])
         c._dependency_link = dependency
    return coeff

%}
