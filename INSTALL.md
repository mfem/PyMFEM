# Installation Guide

## Basic serial install

For serial version, most users (on Linux and Mac) will be fine using the binary bundled in the default `pip` install:

```shell
pip install mfem
```
The above installation will download and install a *serial* version of `MFEM`.


##  Building from source
To build it from source, clone the repository and use pip command

```shell
git clone git@github:mfem/PyMFEM.git
cd PyMFEM
pip install . --user
```

PyMFEM has many options for installation, when building from source, including:
 - Serial and parallel (MPI) wrappers
 - Using pre-built local dependencies
 - Installing additional dependencies such as
   - `hypre`
   - `gslib`
   - `libceed`
   - `metis`
   
Build script checks the following environmental variables
  - CC : c compiler for parallel build
  - CXX : c++ compiler for serial build
  - MPICC : c compiler for parallel build
  - MPICXX : c++ compiler for parallel build
  - CXX11FLAG : C++11 flag for C++ compiler
  - MPIINC : the location of MPI.h (if this variable is set, the parallle PyMFEM is build with CXX, not MPICXX)

Note that `python setup.py install` is deprecated and will be removed soon in favor of `pip`.
When installing via `pip`, options are specified using the `-C` flag using the syntax `-C"name=value"`; e.g. the `--with-parallel` option is now specified as `-C"with-parallel=Yes`.
The name and value of each option should be written explicitly with a dedicated -C flag.

For example, parallel (MPI) support and GSlib support is built with  `--with-parallel`
and `--with-gslib` flags as follows.

```shell
pip install . -C"with-parallel=Yes" -C"with-gslib=Yes"
```

(Warning) The migration to `pip install . ` is on-going effort and Some of the example commands are not tested yet, which are indicated by using old conversion of "python setup.py install XXXX"


## Commonly used flags

| Flag | Description |
|------|-------------|
| `--with-parallel` | Install both serial and parallel versions of `MFEM` and the wrapper<br>(note: this option turns on building `metis` and `hypre`) |
| `--mfem-branch=<reference>` | Download/install MFEM using a specific reference (`git` `branch`, `hash`, or `tag`) |
| `--user` | Install in user's site-package |

## Advanced options

### Suitesparse
`--with-suitesparse` : build MFEM with `suitesparse`. `suitesparse` needs to be installed separately.
Point to the location of `suitesparse` using the flag `--suitesparse-prefix=<location>`

Note: this option turns on building `metis` in serial

### CUDA
`--with-cuda` : build MFEM with CUDA. Hypre cuda build is also supported using
`--with-cuda-hypre`. `--cuda-arch` can be used to specify cuda compute capablility.
(See table in https://en.wikipedia.org/wiki/CUDA#Supported_GPUs)

CUDA needs to be installed separately and nvcc must be found in PATH ([Example](https://github.com/mfem/PyMFEM/blob/e1466a6a/.github/workflows/build-and-test-callable.yml#L111-L122)).

(examples)
```shell
python setup.py install --with-cuda

python setup.py install --with-cuda --with-cuda-hypre

python setup.py install --with-cuda --with-cuda-hypre --cuda-arch=80 (A100)

python setup.py install --with-cuda --with-cuda-hypre --cuda-arch=75 (Turing)
```

### gslib
`--with-gslib` : build MFEM with [GSlib](https://github.com/Nek5000/gslib)

Note: this option downloads and builds GSlib

### libCEED
`--with-libceed` : build MFEM with [libCEED](https://github.com/CEED/libCEED)

Note: this option downloads and builds libCEED

### Specifying compilers
| Flag | Description |
|------|--------|
| `--CC` | c compiler |
| `--CXX` | c++ compiler |
| `--MPICC` | mpic compiler |
| `--MPICXX` | mpic++ compiler |

(example)
Using Intel compiler
```shell
python setup.py install --with-parallel --CC=icc, --CXX=icpc, --MPICC=mpiicc, --MPICXX=mpiicpc
```

### Building MFEM with specific version
By default, setup.py build MFEM with specific SHA (which is usually the released latest version).
In order to use the latest MFEM in Github. One can specify the branch name or SHA using mfem-branch
option.

`--mfem-branch = <branch name or SHA>`

(example)
```shell
pip install . -C"mfem-branch=master"
```

### Using MFEM build externally.
These options are used to link PyMFEM wrapper with existing MFEM library. We need `--mfem-source`
and `--mfem-prefix`

| Flag                       | Description                                                       |
|----------------------------|-------------------------------------------------------------------|
| `--mfem-source <location>` | The location of MFEM source used to build MFEM |
| `--mfem-prefix <location>` | The location of the MFEM library. `libmfem.so` needs to be found in `<location>/lib` |
| `--mfems-prefix <location>`| (optional) Specify serial MFEM location separately |
| `--mfemp-prefix <location>`| (optional) Specify parallel MFEM location separately |


### Blas and Lapack
--with-lapack : build MFEM with lapack

`<location>` is used for CMAKE call to buid MFEM
`--blas-libraries=<location>`
`--lapack-libraries=<location>`

### Options for development and testing
| Flag | Description |
|------|--------|
| `--swig` | run swig only |
| `--skip-swig` | build without running swig` |
| `--skip-ext` | skip building external libraries.|
| `--ext-only` | build exteranl libraries and exit.|

During the development, often we update depenencies (such as MFEM) and edit `*.i` file.
These options allows for building PyMFEM in a step-by-step manner, without rebuilding
the entire wrapper everytime.

First clone the repository or clean everything.

```shell
git clone git@github.com:mfem/PyMFEM.git;cd PyMFEM
python setup.py clean --all
```

Then, build externals alone using the ext-only option
```shell
pip install . -C"ext-only=Yes" --verbose
pip install . -C"with-parallel=Yes" -C"ext-only=Yes" --verbose
```

Then, generate swig wrappers, using the swig option, together with skip-ext, so that
external libraies are not rebuild.
```shell
pip install . -C"skip-ext=Yes"  -C"swig=Yes" --verbose
pip install . -C"with-parallel=Yes" -C"skip-ext=Yes"  -C"swig=Yes" --verbose
```

If you are not happy with the wrapper (`*.cxx` and `*.py`), you edit `*.i` and redo
the same. When you are happy, build the wrapper with skip-swig and skip-ext.

```shell
pip install . -C"skip-ext=Yes"  -C"skip-swig=Yes" --verbose
pip install . -C"with-parallel=Yes" -C"skip-ext=Yes"  -C"skip-swig=Yes" --verbose
```

### Other options
`--unverifiedSSL` :
   This addresses error relating SSL certificate. Typical error message is
   `<urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:xxx)>`
