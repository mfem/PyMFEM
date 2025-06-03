from setuptools import build_meta as _orig
from setuptools.build_meta import *

def get_requires_for_build_wheel(config_settings=None):
    return _orig.get_requires_for_build_wheel(config_settings)


def get_requires_for_build_sdist(config_settings=None):
    return _orig.get_requires_for_build_sdist(config_settings)


def build_wheel(*args, **kwargs):
    return _orig.build_wheel(*args, **kwargs)
