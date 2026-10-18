"""Compare parallel DG advection with the procedural Python ex9p reference.

The reference has fixed module-level options. The fixture changes only these
options in memory, leaving its coefficients, assembly and evolution unchanged.
"""

import importlib.util
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
DATA = EXAMPLES.parent / "data"


def worker(kind, problem, solver):
    from mpi4py import MPI
    problem, solver = int(problem), int(solver)
    if kind == "reference":
        path = EXAMPLES / "ex9p.py"
        source = path.read_text()
        options = dict(problem=problem, ser_ref_levels=1, par_ref_levels=0,
                       order=2, ode_solver_type=solver, t_final=0.05, dt=0.01)
        for name, value in options.items():
            source, count = re.subn(rf"^{name} = .*?$", f"{name} = {value!r}",
                                    source, flags=re.MULTILINE)
            assert count == 1
        # Both path alternatives use this filename in the reference.
        source = source.replace("periodic-hexagon.mesh", "periodic-square.mesh")
        namespace = {"__file__": str(path), "__name__": "ex9p_reference"}
        exec(compile(source, str(path), "exec"), namespace)
        result = dict(solution=namespace["U"].GetDataArray().copy(),
                      steps=namespace["ti"], time=namespace["t"])
    else:
        path = EXAMPLES / "bridge" / "ex9p.py"
        spec = importlib.util.spec_from_file_location("bridge_ex9p", path)
        example = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(example)

        def capture(frame, event, arg):
            if event == "return" and frame.f_code is example.run.__code__:
                for name in ("velocity", "inflow", "u0", "adv"):
                    callbacks = type(frame.f_locals[name])._numba_compiled_overrides
                    assert callbacks and all(callback.nopython_signatures
                                             for callback in callbacks.values())

        sys.setprofile(capture)
        try:
            result = example.run(meshfile=str(DATA / "periodic-square.mesh"),
                                 problem=problem, ode_solver_type=solver,
                                 ser_ref_levels=1, par_ref_levels=0, order=2,
                                 t_final=0.05, dt=0.01)
        finally:
            sys.setprofile(None)
        for options in ({"pa": True}, {"device": "cuda"}):
            try:
                example.run(**options)
            except ValueError as error:
                assert "CPU full assembly only" in str(error)
            else:
                raise AssertionError("unsupported mode accepted")
    np.savez(f"result.{MPI.COMM_WORLD.rank}.npz", **result)


def run_test():
    for problem, solver in [(p, 4) for p in range(4)] + [(0, s) for s in (1, 2, 6)]:
        results = []
        for kind in ("reference", "bridge"):
            with tempfile.TemporaryDirectory(prefix="ex9p-acceptance-") as directory:
                subprocess.run(["mpiexec", "-n", "2", sys.executable, "-I",
                                str(Path(__file__).resolve()), "--worker", kind,
                                str(problem), str(solver)], cwd=directory,
                               check=True, timeout=180)
                ranks = []
                for rank in range(2):
                    with np.load(Path(directory) / f"result.{rank}.npz") as result:
                        ranks.append({name: result[name].copy() for name in result.files})
                results.append(ranks)
        for rank, (expected, actual) in enumerate(zip(*results)):
            assert int(actual["steps"]) == int(expected["steps"]) == 5
            np.testing.assert_allclose(actual["time"], expected["time"], rtol=0, atol=1e-14)
            np.testing.assert_allclose(actual["solution"], expected["solution"],
                                       rtol=1e-10, atol=1e-11)
            error = np.max(np.abs(actual["solution"] - expected["solution"]))
            print(f"PASS problem={problem} solver={solver} rank={rank}: "
                  f"max error={error:.3e}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(*sys.argv[2:])
    else:
        run_test()
