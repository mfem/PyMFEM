# ----------------------------------------------------------------------------------------
# Global build constant parameters
# ----------------------------------------------------------------------------------------
from sys import platform
import os
from shutil import which as find_command
from collections import namedtuple

__all__ = ["swig_command", "rootdir", "extdir",
           "REPOS", "dylibext", "osx_sysroot"]

# ----------------------------------------------------------------------------------------
#  package directory
# ----------------------------------------------------------------------------------------
rootdir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
extdir = os.path.join(rootdir, 'external')
if not os.path.exists(extdir):
    os.mkdir(os.path.join(rootdir, 'external'))

# ----------------------------------------------------------------------------------------
# Platform dependency
# ----------------------------------------------------------------------------------------

osx_sysroot = ''
dylibext = '.so'

if platform == "linux" or platform == "linux2":
    dylibext = '.so'

elif platform == "darwin":
    # OS X
    dylibext = '.dylib'
    import sysconfig
    for i, x in enumerate(sysconfig.get_config_vars()['CFLAGS'].split()):
        if x == '-isysroot':
            osx_sysroot = sysconfig.get_config_vars()['CFLAGS'].split()[i+1]
            break

elif platform == "win32":
    # Windows...
    assert False, "Windows is not supported yet. Contribution is welcome"

# ----------------------------------------------------------------------------------------
# SWIG
# ----------------------------------------------------------------------------------------

swig_command = (find_command('swig') if os.getenv("SWIG") is None
                else os.getenv("SWIG"))
if swig_command is None:
    assert False, "SWIG is not installed (hint: pip install swig)"

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------

release = namedtuple('Release', ['version', 'hash', 'tarball'])
REPOS = dict(
    mfem=dict(
        url="https://github.com/mfem/mfem.git",
        # version, hash, tarball
        releases=[
            release("4.7", "dc9128ef596e84daf1138aa3046b826bba9d259f", None),
            release("4.8", "a01719101027383954b69af1777dc828bf795d62", None),
        ]
    ),
    metis=dict(
        url="https://github.com/KarypisLab/METIS",
        releases=[
            release("5.1.0", "94c03a6e2d1860128c2d0675cbbb86ad4f261256",
                    "https://github.com/mfem/tpls/raw/gh-pages/metis-5.1.0.tar.gz"),
        ]
    ),
    gklib=dict(
        url="https://github.com/KarypisLab/GKlib",
        releases=[
            release("5.1.1", "a7f8172703cf6e999dd0710eb279bba513da4fec",
                    "https://github.com/KarypisLab/GKlib/archive/refs/tags/METIS-v5.1.1-DistDGL-0.5.tar.gz"),
        ]
    ),
    libceed=dict(
        url="https://github.com/CEED/libCEED.git",
        releases=[
            release(
                "0.12.0", None, "https://github.com/CEED/libCEED/archive/refs/tags/v0.12.0.tar.gz"),
        ]
    ),
    hypre=dict(
        url=None,
        releases=[
            release(
                "2.28.0", None, "https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz"),
            release(
                "2.32.0", None, "https://github.com/hypre-space/hypre/archive/v2.32.0.tar.gz"),
        ]
    ),
    gslib=dict(
        url=None,
        releases=[
            release(
                "1.0.8", None, "https://github.com/Nek5000/gslib/archive/refs/tags/v1.0.8.tar.gz"),
        ]
    ),
)
