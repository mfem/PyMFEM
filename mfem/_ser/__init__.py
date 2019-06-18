import mfem

if mfem.mfem_mode is None:
   mfem.mfem_mode = 'serial'

if mfem.mfem_mode == 'parallel':
   raise ImportError("MFEM parallel mode is already loaded")
debug_print = mfem.debug_print


