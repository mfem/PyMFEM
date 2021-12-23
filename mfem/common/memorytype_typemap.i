//
//   this typemap is used to replace MemoryType input
//   to mfem.MemoryType in Pythonl, which is a IntEnum class
//
// usage ENUM_TO_MEMORYTYPE(mfem::MemoryType mt)
%define ENUM_TO_MEMORYTYPE(P1)
  %typemap(in) ( ## P1 ##) {
  PyObject* k = PyObject_GetAttrString($input, "value");
  int i = (int)PyLong_AsLong(k);
  $1 = static_cast< mfem::MemoryType >(i);
}
  
%typemap(typecheck) ( ## P1 ## ){
  $1 = 0;
  PyObject* module = PyImport_ImportModule("enum");
  if (!module){
      $1 = 0;
  } else {     
      PyObject* cls = PyObject_GetAttrString(module, "IntEnum");
      if (!cls){
          $1 = 0;            
      } else {
         int check = PyObject_IsInstance($input, cls);
	 if (check) {
	   $1 = 1;
	 }
         Py_DECREF(cls);	 
      }
      Py_DECREF(module);
  }
}
%enddef
