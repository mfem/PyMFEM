# Installation Guide

## Basic install

Most users will be fine using the binary bundled in the default `pip` install:

```shell
pip install mfem
```
The above installation will download and install a *serial* version of `MFEM`.

##  Building from source
PyMFEM has many options for installation, when building from source, including:
 - Serial and parallel (MPI) wrappers
 - Using pre-built local dependencies
 - Installing additional dependencies such as
   - `hypre`
   - `gslib`
   - `libceed`
   - `metis`

Most of the options for PyMFEM can be used directly when installing via `python setup.py install`, e.g.
```shell
git clone git@github:mfem/PyMFEM.git
cd PyMFEM
python setup.py install --user
```
For example, parallel (MPI) support is built with  the `--with-parallel` flag:
```shell
python setup.py install --with-parallel
```

Note: this option turns on building `metis` and `Hypre`

## Commonly used flags

| Flag | Description |
|------|-------------|
| `--with-parallel` | Install both serial and parallel versions of `MFEM` and the wrapper<br>(note: this option turns on building `metis` and `hypre`) |
| `--mfem-branch=<reference>` | Download/install MFEM using a specific reference (`git` `branch`, `hash`, or `tag`) |
| `--mfem-source=<location>` | Specify a local version of MFEM |
| `--mfem-prefix <location>` | The location of the MFEM library. `libmfem.so` needs to be found in `<location>/lib` |
| `--mfems-prefix <location>`| (optional) Specify serial MFEM location separately |
| `--mfemp-prefix <location>`| (optional) Specify parallel MFEM location separately |

In order to see the full list of options, use

```shell
python setup.py install --help
```

## Advanced options

### Suitesparse
`--with-suitesparse` : build MFEM with `suitesparse`. `suitesparse` needs to be installed separately.
Point to the location of `suitesparse` using the flag `--suitesparse-prefix=<location>`

Note: this option turns on building `metis` in serial

### CUDA
`--with-cuda` : build MFEM with CUDA. Hypre cuda build is also supported using
`--with-cuda-hypre`. `--cuda-arch` can be used to specify cuda compute capablility.
(See table in https://en.wikipedia.org/wiki/CUDA#Supported_GPUs)

CUDA needs to be installed separately and nvcc must be found in PATH.

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
python setup.py install --mfem-branch=master
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

### Options fordDevelopment and testing
`--swig` : run swig only
`--skip-swig` : build without running swig
`--skip-ext` : skip building external libraries.
`--ext-only` : build exteranl libraries and exit.

During the development, often we update depenencies (such as MFEM) and edit `*.i` file.

First clean everything.

```shell
python setup.py clean --all
```

Then, build externals alone
```shell
python setup.py install --with-parallel --ext-only --mfem-branch=master
```

Then, genrate swig wrappers.
```shell
python setup.py install --with-parallel --swig --mfem-branch=master
```

If you are not happy with the wrapper (`*.cxx` and `*.py`), you edit `*.i` and redo
the same. When you are happy, build the wrapper. `--swig` does not clean the
existing wrapper. So, it will only update wrapper for updated `*.i`

When building a wrapper, you can use `--skip-ext` option. By default, it will re-run
swig to generate entire wrapper codes.
```shell
python setup.py install --with-parallel --skip-ext --mfem-branch=master
```

If you are sure, you could use `--skip-swig` option, so that it compiles the wrapper
codes without re-generating it.
```shell
python setup.py install --with-parallel --skip-ext --skip-swig --mfem-branch=master
```

### Other options
`--unverifiedSSL` :
   This addresses error relating SSL certificate. Typical error message is
   `<urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:xxx)>`


