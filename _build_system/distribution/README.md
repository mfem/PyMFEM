Build and validate the serial and MPICH wheels locally using the same scripts
as CI. Run from the repository root with Python 3.11–3.14.

On Linux, install a C/C++ compiler, `readelf`, and zlib development headers.
For the CI platform, use `quay.io/pypa/manylinux_2_28_x86_64`.
On macOS, install Xcode command-line tools with `xcode-select --install`.

Create and activate a clean build environment, then install the tools and MPI:

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
bash _build_system/distribution/build_serial_wheel.sh \
    --work-dir /tmp/pymfem-serial-candidate

bash _build_system/distribution/build_mpich_wheel.sh \
    --work-dir /tmp/pymfem-mpich-candidate \
    --base-wheel /tmp/pymfem-serial-candidate/mfem-*.whl
```

Use new work directories outside the checkout. Each script copies the checkout,
builds and validates there, and places the final wheel and `build.log` directly
in its work directory. Builds default to two workers; set `PYMFEM_BUILD_NCPU`
to change this.
macOS builds use the native architecture and target macOS 14 by default
(`MACOSX_DEPLOYMENT_TARGET` overrides it). Apple Silicon is the prepared CI
target; native macOS build and runtime validation remain pending.

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

The validator creates the runtime venv, installs dependencies from PyPI, and
checks serial/parallel solves, GSLIB, MPI, and uninstall/reinstall behavior.
Results are saved in `/tmp/pymfem-downloaded-validation/validation.json`.

To test upgrading from a previously validated candidate, keep the two versions
in separate wheel directories and run:

```bash
python _build_system/distribution/validate_mpich_upgrade.py \
    /tmp/pymfem-previous-wheels /tmp/pymfem-downloaded \
    /tmp/pymfem-upgrade-validation
```

This installs and upgrades through `mfem[mpich]`, runs the solves, and checks
the diagnostic for mismatched versions before restoring the matching pair.
