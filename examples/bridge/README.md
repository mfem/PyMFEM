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

## Standalone example 9

`ex9.py` implements the serial DG advection problem without importing the
ordinary ex9 driver. Velocity, initial value and zero inflow are native NSB
coefficient directors. A native TimeDependentOperator callback computes
`M^{-1}(K u + b)` through SparseMatrix and CGSolver ordinary bridges.
Matrices, forms, the mass smoother and solver are configured and retained in
Python; reusable callback scratch assumes non-reentrant execution.

Requires an ABI-7 bridge-enabled PyMFEM build and `nsb_rt>=0.15.0`:

```shell
python examples/bridge/ex9.py -r 0 -o 1 -p 0 -s 4 -tf 0.02 -dt 0.01
```

The example supports CPU full assembly and explicit solvers 1, 2, 3, 4 and 6.
Partial assembly, other devices and implicit schemes are rejected. Problems
2 and 3 require dimension at least two. Bounds/problem checks precede native
callback compilation. Files and visualization options follow ordinary ex9.
Output reports the number of unknowns and time-step progress, as in ordinary
ex9. The distributed example contains no performance instrumentation.

## Standalone example 19

`ex19.py` implements serial incompressible nonlinear elasticity without
importing the ordinary example. NSB compiles the deformation/reference
coefficients, the rubber operator's `Mult` and `GetGradient`, and the block
preconditioner's `SetOperator` and `Mult`. Solver setup and residual monitors
remain in Python. It requires the same ABI-7 runtime as ex9 and currently
uses CPU execution.

```shell
python examples/bridge/ex19.py -m beam-quad.mesh -r 0 -o 2
```

The nonlinear form owns the Jacobian. Each gradient rebuild invalidates
previous borrowed references; Newton updates the preconditioner before its
next application. The example retains the form, pressure mass matrix,
solvers, and smoothers throughout the solve. Bounded vector views borrow
callback storage locally; reusable scratch assumes non-reentrant execution.
Visualization is opt-in with `-vis`. Output files follow ordinary ex19.
