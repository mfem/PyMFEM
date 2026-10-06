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

`test_coefficient_bridges.py` exercises direct `Coefficient`,
`VectorCoefficient`, and `MatrixCoefficient` NSB directors through MFEM
projection and diffusion assembly. Example-4 comparisons and timing commands
are recorded in the shared integration report
[NSBR_05](../../../dev_guide/nsb_pymfem_integration/NSBR_05_PyMFEM_ex4_ex38_application.md).

`test_data_members.py` verifies automatically bridged IntegrationPoint fields,
ElementTransformation element numbers, and const integration-point reads through
MFEM coefficient projection. It also rejects writes to the const callback input.
