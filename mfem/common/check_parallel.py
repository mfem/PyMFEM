"""
Check matching versions of separately installed serial and parallel wheels.
"""

from importlib.metadata import PackageNotFoundError, version


def check_parallel_version():
    """
    Reject mismatched binary distributions before loading native extensions.
    """
    try:
        addon = version('mfem-mpich')
    except PackageNotFoundError:
        return  # A conventional source build has no separate add-on.
    try:
        base = version('mfem')
    except PackageNotFoundError:
        base = 'missing'
    if base != addon:
        raise ImportError(
            f'Incompatible PyMFEM wheels: mfem={base}, mfem-mpich={addon}. '
            'Install matching versions with pip install --upgrade "mfem[mpich]".')
