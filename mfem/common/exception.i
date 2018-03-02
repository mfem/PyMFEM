%exception {
  try {
    $action
  }
  %#ifdef  MFEM_USE_EXCEPTIONS
  catch (mfem::ErrorException &_e) {
         std::string s("PyMFEM error (mfem::ErrorException): "), s2(_e.what());
         s = s + s2;    
         SWIG_exception(SWIG_RuntimeError, s.c_str());
  }
  %#endif

  catch (...) {
        SWIG_exception(SWIG_RuntimeError, "unknown exception");
  }	 
}
