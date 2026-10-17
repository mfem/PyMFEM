"""Compare native serial ex19 directors with ordinary Python Newton solves."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import numpy as np


EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
DATA = EXAMPLES.parent / "data"


def worker(kind, mesh):
    path = EXAMPLES / "ex19.py" if kind == "reference" else EXAMPLES / "bridge" / "ex19.py"
    spec = importlib.util.spec_from_file_location("ex19_under_test", path)
    example = importlib.util.module_from_spec(spec)
    # Ordinary ex19 parses arguments at import time.
    sys.argv = [str(path)]
    spec.loader.exec_module(example)
    captured = {}

    def capture(frame, event, arg):
        if event != "return":
            return
        if frame.f_code.co_name == "Solve" and frame.f_globals is example.__dict__:
            operator = frame.f_locals["self"]
            solver = operator.newton_solver
            captured["solution"] = frame.f_locals["xp"].GetDataArray().copy()
            captured["iterations"] = solver.GetNumIterations()
            captured["initial_norm"] = solver.GetInitialNorm()
            captured["final_norm"] = solver.GetFinalNorm()
            captured["converged"] = solver.GetConverged()
            if kind == "bridge":
                for director in (operator, operator.j_prec):
                    callbacks = type(director)._numba_compiled_overrides
                    assert callbacks and all(callback.nopython_signatures
                                             for callback in callbacks.values())
        if kind == "bridge" and frame.f_code is example.ex19_main.__code__:
            for name in ("deform", "refconfig"):
                if name not in frame.f_locals:
                    continue
                callbacks = type(frame.f_locals[name])._numba_compiled_overrides
                assert callbacks and all(callback.nopython_signatures
                                         for callback in callbacks.values())

    options = SimpleNamespace(mesh=str(DATA / mesh), refine=0, order=2,
                              relative_tolerance=1e-4, absolute_tolerance=1e-6,
                              newton_iterations=500, shear_modulus=1.0,
                              visualization=False)
    sys.setprofile(capture)
    try:
        example.ex19_main(options)
    finally:
        sys.setprofile(None)
    assert captured["converged"]
    assert captured["final_norm"] <= max(1e-6, 1e-4 * captured["initial_norm"])
    np.savez("result.npz", **captured)


def run_test():
    for mesh in ("beam-quad.mesh", "beam-tri.mesh", "beam-hex.mesh"):
        results = []
        for kind in ("reference", "bridge"):
            with tempfile.TemporaryDirectory(prefix="ex19-acceptance-") as directory:
                subprocess.run([sys.executable, "-I", str(Path(__file__).resolve()),
                                "--worker", kind, mesh], cwd=directory, check=True)
                with np.load(Path(directory) / "result.npz") as result:
                    results.append({name: result[name].copy() for name in result.files})
        expected, actual = results
        assert int(actual["iterations"]) == int(expected["iterations"])
        np.testing.assert_allclose(actual["solution"], expected["solution"],
                                   rtol=1e-8, atol=1e-10)
        error = np.max(np.abs(actual["solution"] - expected["solution"]))
        print(f"PASS {mesh}: iterations={int(actual['iterations'])}, "
              f"final norm={float(actual['final_norm']):.3e}, max error={error:.3e}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(*sys.argv[2:])
    else:
        run_test()
