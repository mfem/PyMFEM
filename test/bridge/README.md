# PyMFEM bridge tests

Build and install PyMFEM with the serial Numba bridge enabled, then run:

```bash
python test/test_bridge.py
```

`test_bridge.py` discovers and runs each `test_*.py` file in this directory.
These tests are deliberately excluded from `test/run_tests.py`: an ordinary
PyMFEM build does not include the generated bridge registration module.

Each bridge test should verify numerical results and that its call path compiles
in `numba.njit`. Keep performance measurements in `examples/bridge/`, where
their timing output is useful without making a test flaky.

`test_ex9.py` compares five-step serial trajectories with ordinary ex9 in
separate installed-package processes. It covers all four 2D problems, the 1D
problems, and all supported explicit solvers, checks compiled callback
signatures, and verifies rejection of partial assembly and device execution.

`test_solvers.py` checks ordinary solver registration, absence of factories
for abstract solver bases, inherited configuration calls on six concrete
solvers, and construction of FGMRESSolver inside compiled code.

`test_ex19.py` compares ordinary and bridged Newton solves on quadrilateral,
triangular, and hexahedral beam meshes, checking convergence, iteration counts,
solution fields, and compiled callback signatures.

`test_coefficient_bridges.py` exercises direct `Coefficient`,
`VectorCoefficient`, and `MatrixCoefficient` NSB directors through MFEM
projection and diffusion assembly. Example-4 comparisons and timing commands
are recorded in the shared integration report
[NSBR_05](../../../dev_guide/nsb_pymfem_integration/NSBR_05_PyMFEM_ex4_ex38_application.md).

`test_data_members.py` verifies automatically bridged IntegrationPoint fields,
ElementTransformation element numbers, and const integration-point reads through
MFEM coefficient projection. It also rejects writes to the const callback input.

`test_linearform_integrators.py` exercises direct `LinearFormIntegrator` NSB
directors through `LinearForm.Assemble()`, including native ownership transfer,
facade closure before assembly, surface and subdomain shape paths, and comparison
with Python callbacks and MFEM's native `DomainLFIntegrator`. Its callbacks and
the direct coefficient tests also exercise explicit callback-local
`stateaccess` declarations through real MFEM projection and assembly calls.


## Parallel foundation (NSBP_13 Stage 4)

Build with both `-C"with-parallel=Yes"` and
`-C"with-numba-swig-bridge=Yes"`. Parallel artifacts require NSB runtime
0.15.1 or later. Serial and parallel imports must use separate processes.
MPI communicators and native objects stay local to each rank; applications
keep borrowed MPI resources valid until native use and teardown have ended.

Run each check with a timeout and at least two ranks:

```sh
mpiexec -n 2 python test/bridge/test_parallel_foundation.py
mpiexec -n 2 python test/bridge/test_ex18p.py
mpiexec -n 2 python examples/bridge/ex4p.py --refine 0 --compare
mpiexec -n 2 python examples/bridge/ex18p.py -rs 0 -rp 0 -o 1 -tf 0.02
```

`test_parallel_foundation.py` checks ordinary ownership handoff, Hypre class
coverage, nested directors, rank-local callback errors and recovery, teardown,
and native `GetComm()` passed to a compiled CGSolver constructor. The ex18p
fixture compares the rank-local RHS and characteristic speed with the Python
reference. The ex4p comparison checks local RHS, projection, solution and global
L2 error. These are CPU checks; MPI collectives remain in the Python driver or
native MFEM methods. No ex38p translation is included.

Developer regeneration uses `generate_nsb_artifacts(serial=False)` from
`_build_system/build_pymfem.py`, with `build_globals.hypre_prefix` set to the
configured Hypre installation. It reads MPI-enabled MFEM headers from
`external/mfem/cmbuild_par`. The helper supplies mpi4py and Hypre SWIG includes
without storing workstation paths in `bridge.toml` or the manifest. Pass
`check=True` for the parallel artifact freshness check. Consumer builds use
committed artifacts and need only the runtime, MPI/mpi4py and native build tools.

## Parallel advection and nonlinear elasticity (Stage 5)

The acceptance drivers launch their own two-rank MPI workers. Run them once
from the installed consumer environment (do not wrap these drivers in mpiexec):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I test/bridge/test_ex9p.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I test/bridge/test_ex19p.py
```

ex9p compares four advection problems and explicit Euler/RK2/RK4/RK6 on a
periodic square mesh. ex19p compares Newton solves on beam quad, tri and hex
meshes, including repeated Jacobian/preconditioner updates and rank-local
solutions. Both verify compiled callback signatures against ordinary Python
references in separate processes. The examples can also be run directly:

```sh
mpiexec -n 2 python examples/bridge/ex9p.py -rs 0 -rp 0 -o 2 -tf 0.05
mpiexec -n 2 python examples/bridge/ex19p.py -m beam-quad.mesh
```

ex9p supports CPU full assembly; partial assembly and GPU modes are rejected.
ex19p retains its pressure mass OperatorHandle as the single matrix owner.
