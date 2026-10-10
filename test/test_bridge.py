"""Run the serial numba-swig-bridge integration tests.

This suite is intentionally separate from the ordinary PyMFEM tests because it
requires a PyMFEM installation built with ``with-numba-swig-bridge=Yes``.
"""
from __future__ import print_function

from pathlib import Path
import subprocess
import sys


TEST_DIR = Path(__file__).resolve().parent
BRIDGE_DIR = TEST_DIR / "bridge"


def run_test():
    tests = sorted(BRIDGE_DIR.glob("test_*.py"))
    if not tests:
        raise RuntimeError("no bridge tests found")

    failed = []
    for test in tests:
        print("#### running bridge test: " + test.name)
        result = subprocess.run([sys.executable, str(test)], check=False)
        if result.returncode:
            failed.append(test.name)

    if failed:
        raise SystemExit("bridge tests failed: " + ", ".join(failed))


if __name__ == "__main__":
    if "-p" in sys.argv[1:]:
        print("bridge tests currently cover mfem.ser only; skipping parallel suite")
    else:
        run_test()
