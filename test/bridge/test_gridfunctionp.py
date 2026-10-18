"""Two-rank parallel variant of the shared native bridge fixture."""

from pathlib import Path
import runpy
import sys

sys.argv.append("-p")
runpy.run_path(
    str(Path(__file__).with_name("test_gridfunction.py")), run_name="__main__"
)
