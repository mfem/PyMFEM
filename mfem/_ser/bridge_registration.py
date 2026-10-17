"""Install the serial numba-swig-bridge registrations for this PyMFEM build."""
try:
    from . import _nsb_registration
except ImportError as error:
    if error.name != "mfem._ser._nsb_registration":
        raise
    raise ImportError(
        "This PyMFEM installation was built without numba-swig-bridge. "
        'Reinstall with -C"with-numba-swig-bridge=Yes".'
    ) from error
