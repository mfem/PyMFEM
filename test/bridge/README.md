# PyMFEM bridge tests

## Test All cases
Run `python ../test_bridge.py`

## Test cases
### test_coefficient_bridges.py

Checks coefficient directors through MFEM projection and diffusion assembly.

### test_data_members.py

Checks public C++ data members and rejects writes through const callback
arguments.

### test_densemat.py

Checks DenseMatrix methods called from Numba-compiled code.

### test_densetensor.py

Checks DenseTensor methods and zero-copy data views.

### test_ex18.py

Compares the serial ex18 compiled RHS with the Python reference.

### test_ex18p.py

Compares the parallel ex18 RHS and characteristic speed with the Python
reference. Run on two MPI ranks.

### test_ex19.py

Compares serial Newton solves, solution fields, convergence and callback
compilation on three beam meshes.

### test_ex19p.py

Compares parallel Newton solves with the Python reference on three beam meshes.
Launches its own two-rank workers.

### test_ex9.py

Compares short serial advection runs across problems and explicit solvers;
checks compiled callbacks and unsupported modes.

### test_ex9p.py

Compares parallel advection problems and explicit solvers with the Python
reference. Launches its own two-rank workers.

### test_fespace.py

Checks finite-element-space methods and callback-scoped borrowed returns.

### test_hyperbolic.py

Checks FluxFunction and HyperbolicFormIntegrator methods called from compiled
code.

### test_intarray.py

Checks the ordinary bridge for Array<int>.

### test_linearform_integrators.py

Checks integrator directors through linear-form assembly, including native
ownership and comparisons with Python and native integrators.

### test_nonlinearform.py

Checks NonlinearForm methods called from compiled code.

### test_nsb_state_build.py

Checks compiler optimization flags for native state sources, including
concurrent builds. This is a source-level build test.

### test_parallel_foundation.py

Checks Hypre bridges, nested directors, communicator handling, ownership and
callback error recovery. Run on two MPI ranks.

### test_solvers.py

Checks concrete and abstract solver registration, inherited methods and
construction in compiled code.

### test_time_dependent_operator.py

Checks native C++ dispatch to a compiled void-return director override.

### test_vector.py

Checks Vector methods and zero-copy writes from compiled code.

### Native operators and GridFunction

`test_native_operators.py` checks native operators, merged vector/matrix Add,
both vector Subtract overloads in Python and compiled code,
allocating helpers and BlockMatrix Mult, and rejection of deferred/reflected
forms. `test_gridfunction.py` checks GridFunction methods, Vector inheritance,
compiled construction/ownership transfer, and direct state-field compound
assignment. Run each script with `-p` for the
parallel bindings; the GridFunction script also exercises ParGridFunction and
must be run under MPI for parallel coverage. Keep serial and parallel runs in
separate Python processes.
