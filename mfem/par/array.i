%module array
%rename(Equal) mfem::Array <class T>::operator=;
%{
#include "general/array.hpp"
#include <iostream>
#include "numpy/arrayobject.h"    
%}

%init %{
import_array();
%}

%typemap(in) (int *_data, int asize) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $2 = PyList_Size($input);
  $1 = (int *) malloc(($2)*sizeof(int));
  for (i = 0; i < $2; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyInt_Check(s)) {
        free($1);
        PyErr_SetString(PyExc_ValueError, "List items must be integer");
        return NULL;
    }
    $1[i] = (int)PyInt_AsLong(s);
  }
}
%typemap(typecheck) (int *_data, int asize) {
   $1 = PyList_Check($input) ? 1 : 0;
}

%typemap(newfree) (int *_data,  int asize) {
   if ($1) free($1);
}

//%import "intrules.i"
//%newobject intArray
%include "general/array.hpp"
%extend mfem::Array{
  void __setitem__(int i, const T v) {
    (* self)[i] = v;
    }
  const T & __getitem__(const int i) const{
    return (* self)[i];
  }
  void Assign(const T &a){
     *self = a;
  }   
  void FakeToList(void){}
};
namespace mfem{
%feature("shadow")Array::FakeToList %{
def ToList(self):
    return [self[i] for i in range(self.Size())]
%}
}  
namespace mfem{
%template(intArray) Array<int>;
}
namespace mfem{
%template(doubleArray) Array<double>;
}


 /*
%pythoncode %{
def ToList(A):
    return [A[i] for i in range(A.Size())]
%}


namespace mfem{
class BaseArray
{
protected:
   /// Pointer to data
   void *data;
   /// Size of the array
   int size;
   /// Size of the allocated memory
   int allocsize;
   int inc;
   
   BaseArray() { }
   /// Creates array of asize elements of size elementsize
   BaseArray(int asize, int ainc, int elmentsize);
   /// Free the allocated memory
   ~BaseArray();
   void GrowSize(int minsize, int elementsize);
};
%rename(intArray) Array<int>; 

class Array<int> : public BaseArray
{
 public:  
     /// Creates array of asize elements
   explicit inline Array(int asize = 0, int ainc = 0)
      : BaseArray(asize, ainc, sizeof (int)) { }

   inline Array(int *_data, int asize, int ainc = 0)
     { data = _data; size = asize; allocsize = -asize; inc = ainc; }
   /// Destructor
   inline ~Array() { }
   
   inline int Size() const { return size; }
   /// Prints array to stream with width elements per row
   void Print(std::ostream &out = std::cout, int width = 4);
};

}
*/
