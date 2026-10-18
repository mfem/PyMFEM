"""Run numba-swig-bridge integration tests (all by default, -p for parallel).

This suite is intentionally separate from the ordinary PyMFEM tests because it
requires a PyMFEM installation built with ``with-numba-swig-bridge=Yes``.
"""
from __future__ import print_function

from pathlib import Path
import argparse
import shutil
import subprocess
import sys


TEST_DIR = Path(__file__).resolve().parent
BRIDGE_DIR = TEST_DIR / "bridge"
MPI_TESTS = {"test_parallel_foundation.py", "test_ex18p.py",
             "test_gridfunctionp.py", "test_native_operatorsp.py"}
PARALLEL_TESTS = MPI_TESTS | {"test_ex9p.py", "test_ex19p.py"}


def run_test(mode="all"):
    tests = sorted(BRIDGE_DIR.glob("test_*.py"))
    if mode == "parallel":
        tests = [test for test in tests if test.name in PARALLEL_TESTS]
    elif mode == "serial":
        tests = [test for test in tests if test.name not in PARALLEL_TESTS]
    if not tests:
        raise RuntimeError("no bridge tests found")

    launcher = None
    if any(test.name in MPI_TESTS for test in tests):
        launcher = shutil.which("mpiexec")
        if launcher is None:
            raise RuntimeError("MPI launcher unavailable: mpiexec")

    failed = []
    for test in tests:
        print("#### running bridge test: " + test.name, flush=True)
        command = [sys.executable, str(test)]
        if test.name in MPI_TESTS:
            command = [launcher, "-n", "2", *command]
        result = subprocess.run(command, check=False)
        if result.returncode:
            failed.append(test.name)

    if failed:
        raise SystemExit("bridge tests failed: " + ", ".join(failed))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("-p", "--parallel", action="store_true")
    modes.add_argument("--serial", action="store_true")
    args = parser.parse_args()
    mode = "parallel" if args.parallel else "serial" if args.serial else "all"
    run_test(mode)
