from setuptools import build_meta as _orig
from setuptools.build_meta import *

import build_globals as bglb


def _enabled(config_settings, flag):
    """Return whether a PEP 517 boolean build setting is enabled."""
    if not config_settings:
        return False
    value = config_settings.get(flag, "No")
    if isinstance(value, (list, tuple)):
        value = value[-1] if value else "No"
    return str(value).upper() in ("YES", "TRUE", "1")


def get_requires_for_build_wheel(config_settings=None):
    # Do not consume config_settings here: build_wheel passes the same setting
    # on to setup.py, where build_config.py enables the matching feature.
    ret = _orig.get_requires_for_build_wheel(config_settings)
    if _enabled(config_settings, "with-parallel"):
        ret = ret + ['mpi4py']
    if _enabled(config_settings, "with-numba-swig-bridge"):
        ret = ret + ['numba-swig-bridge-rt>=0.14.0']
    return ret


def get_requires_for_build_sdist(config_settings=None):
    return _orig.get_requires_for_build_sdist(config_settings)


def build_wheel(*args, **kwargs):
    bglb.cfs = args[1]
    if bglb.cfs is None:
        bglb.cfs = {}
    else:
        # process_cmd_options consumes its input. Keep the caller-owned PEP 517
        # mapping intact for build frontends that reuse it.
        bglb.cfs = dict(bglb.cfs)
    return _orig.build_wheel(*args, **kwargs)
