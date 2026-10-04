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
