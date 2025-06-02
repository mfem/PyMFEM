from setuptools import build_meta as _orig
from setuptools.build_meta import *
import build_utils

build_config_options = None


def get_requires_for_build_wheel(config_settings=None):
    build_options = {}
    if config_settings is not None:
        build_options = config_settings.pop("build_options", {})

    ret = _orig.get_requires_for_build_wheel(config_settings)
    return ret


def get_requires_for_build_sdist(config_settings=None):
    return _orig.get_requires_for_build_sdist(config_settings)


def build_wheel(*args, **kwargs):
    return _orig.build_wheel(*args, **kwargs)
