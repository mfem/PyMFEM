Build and validate the serial and MPICH wheels locally using the same scripts
as CI. Run from the repository root with Python 3.11–3.14.
Create and activate a clean build environment, then install the tools and MPI:

## TL;DR

```bash
python -m venv /tmp/build_serial
source /tmp/build_serial/bin/activate
pip install --upgrade pip
pip install wheel packaging twine swig
```

Install the repair tool for your platform, then MPI:

```bash
# Linux:
pip install auditwheel patchelf
# macOS:
pip install 'delocate==0.13.0'

pip install 'mpich==5.0.2' 'mpi4py==4.1.2'
```

Set the desired version in `mfem/__init__.py`, then build both wheels:

```bash
version=$(python -c 'import mfem; print(mfem.__version__)')
# create serial wheel. pass workdirectory where serial wheel will be created.
bash _build_system/distribution/build_serial_wheel.sh \
    /tmp/pymfem-serial-candidate --plat PLATFORM

# create parallel wheel. pass workdirectory and serial wheel created abovel.
bash _build_system/distribution/build_mpich_wheel.sh \
    /tmp/pymfem-mpich-candidate \
    /tmp/pymfem-serial-candidate/mfem-*.whl --plat PLATFORM
```

After both scripts succeed, upload manually to TestPyPI using a fresh version:

```bash
twine upload --repository testpypi \
    /tmp/pymfem-serial-candidate/*.whl \
    /tmp/pymfem-mpich-candidate/*.whl
```

Download that version from TestPyPI and validate it in a fresh runtime venv:

```bash
pip download --index-url https://test.pypi.org/simple/ \
    --no-deps --only-binary=:all: --dest /tmp/pymfem-downloaded \
    "mfem==$version" "mfem-mpich==$version"

python _build_system/distribution/validate_mpich_pair.py \
    /tmp/pymfem-downloaded/mfem-*.whl \
    /tmp/pymfem-downloaded/mfem_mpich-*.whl \
    --work-dir /tmp/pymfem-downloaded-validation \
    --examples-dir examples --mesh data/star.mesh
```


To test upgrading from a previously validated candidate, keep the two versions
in separate wheel directories and run:

```bash
python _build_system/distribution/validate_mpich_upgrade.py \
    /tmp/pymfem-previous-wheels /tmp/pymfem-downloaded \
    /tmp/pymfem-upgrade-validation
```

This installs and upgrades through `mfem[mpich]`, runs the solves, and checks
the diagnostic for mismatched versions before restoring the matching pair.

## Extra Note

### Platform selection
On Linux, PLATFORM default is `manylinux_2_28` on hosts with glibc 2.28 or older.
On newer hosts, they use auditwheel's `auto` mode to select a compatible tag
from the built wheel.
In the current release CI, `manylinux_2_28_x86_64` and `quay.io/pypa/manylinux_2_28_x86_64`
is used.
On macOS, install Xcode command-line tools with `xcode-select --install`.

### Build details
Use new work directories outside the checkout. Each script copies the current
sources, including local edits, without build caches, compiled libraries,
generated SWIG wrappers, or downloaded external dependencies. It builds and
validates there, and places the final wheel and `build.log` directly
in its work directory. Builds default to two workers; set `PYMFEM_BUILD_NCPU`
to change this.
macOS builds use the native architecture and target macOS 14 by default
(`MACOSX_DEPLOYMENT_TARGET` overrides it). Apple Silicon is the prepared CI
target; native macOS build and runtime validation remain pending.

### Validation
The validator creates the runtime venv, installs dependencies from PyPI, and
checks serial/parallel solves, GSLIB, MPI, and uninstall/reinstall behavior.
Constraints keep `mfem` and `mfem-mpich` on the supplied wheel files during
extra installation and upgrade validation, even if PyPI has the same versions.
Results are saved in `/tmp/pymfem-downloaded-validation/validation.json`.
