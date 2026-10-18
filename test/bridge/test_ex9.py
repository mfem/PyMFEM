"""Compare installed serial ex9 callbacks with the ordinary Python example."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
DATA = EXAMPLES.parent / "data"


def worker(kind, mesh, problem, solver):
    path = EXAMPLES / "ex9.py" if kind == "reference" else EXAMPLES / "bridge" / "ex9.py"
    spec = importlib.util.spec_from_file_location("ex9_under_test", path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    captured = {}

    # The ordinary example does not return its solution. Capture the final
    # field before its run frame releases the native owners, without changing it.
    def capture(frame, event, arg):
        if event == "return" and frame.f_code is example.run.__code__:
            captured["solution"] = frame.f_locals["u"].GetDataArray().copy()
            captured["steps"] = frame.f_locals["ti"]
            captured["time"] = frame.f_locals["t"]
            if kind == "bridge":
                for name in ("velocity", "inflow", "u0", "adv"):
                    callbacks = type(frame.f_locals[name])._numba_compiled_overrides
                    assert callbacks and all(callback.nopython_signatures
                                             for callback in callbacks.values())

    sys.setprofile(capture)
    try:
        example.run(meshfile=str(DATA / mesh), problem=int(problem),
                    ode_solver_type=int(solver), ref_levels=0, order=2,
                    t_final=0.05, dt=0.01, visualization=False)
    finally:
        sys.setprofile(None)
    if kind == "bridge":
        for options in ({"pa": True}, {"device": "cuda"}):
            try:
                example.run(**options)
            except ValueError as error:
                assert "CPU full assembly only" in str(error)
            else:
                raise AssertionError("unsupported execution mode was accepted")
    np.savez("result.npz", **captured)


def run_test():
    cases = [("periodic-square.mesh", p, 4) for p in range(4)]
    cases += [("periodic-segment.mesh", p, 4) for p in (0, 1)]
    cases += [("periodic-square.mesh", 0, s) for s in (1, 2, 3, 6)]
    for mesh, problem, solver in cases:
        results = []
        for kind in ("reference", "bridge"):
            # Separate processes isolate MFEM's global Device configuration.
            with tempfile.TemporaryDirectory(prefix="ex9-acceptance-") as directory:
                subprocess.run([sys.executable, "-I", str(Path(__file__).resolve()),
                                "--worker", kind, mesh, str(problem), str(solver)],
                               cwd=directory, check=True)
                with np.load(Path(directory) / "result.npz") as result:
                    results.append({name: result[name].copy() for name in result.files})
        expected, actual = results
        assert int(actual["steps"]) == int(expected["steps"]) == 5
        np.testing.assert_allclose(actual["time"], expected["time"], rtol=0, atol=1e-14)
        np.testing.assert_allclose(actual["solution"], expected["solution"],
                                   rtol=1e-10, atol=1e-11)
        error = np.max(np.abs(actual["solution"] - expected["solution"]))
        print(f"PASS {mesh} problem={problem} solver={solver}: max error={error:.3e}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(*sys.argv[2:])
    else:
        run_test()
