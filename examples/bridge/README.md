# Numba bridge examples

`vector.py` is the first `numba-swig-bridge` example for PyMFEM. It imports the
serial Vector registration generated during a bridge-enabled build, then calls
`mfem::Vector::Norml2()` repeatedly from an `@numba.njit` function. The public
`mfem.get_bridge_registration()` installs and returns the enabled serial bridge registration.
Use the public imports below so the local ``mfem`` name continues to refer to
the serial interface:

```python
import mfem.ser as mfem
registration = mfem.get_bridge_registration()
```

Build PyMFEM with the bridge enabled, then run the example:

```shell
python -m pip install ".[bridge]" --no-build-isolation \
  -C"with-numba-swig-bridge=Yes"
python examples/bridge/vector.py
```

The example compiles the Numba function before timing. It then runs the same
`Norml2` loop through the Numba bridge and through PyMFEM's usual Python proxy,
and reports both times and their ratio. Use `--count` to change the number of
calls, for example `python examples/bridge/vector.py --count 10000`.

`ex18.py` applies the same technique to MFEM example 18. It keeps setup and
time integration in Python, while a Numba director performs the per-element
volume contribution using `GetDataArray()` views. It currently supports the
preassembled weak-divergence path and uniform-order spaces:

```shell
python examples/bridge/ex18.py -r 0 -o 1 -tf 0.02
```
