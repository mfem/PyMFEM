"""Check ordinary solver registration and inherited compiled method calls."""

import json
from pathlib import Path

import numpy as np
from numba import njit
import nsb_rt as nsb
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def solve(solver, matrix, rhs, solution):
    # These configuration methods belong to the abstract IterativeSolver base.
    solver.SetRelTol(1e-12)
    solver.SetAbsTol(1e-12)
    solver.SetMaxIter(100)
    solver.SetPrintLevel(0)
    solver.SetOperator(matrix)
    solver.Mult(rhs, solution)
    return solver.GetConverged()


@njit
def make_fgmres():
    return nsb.own_this(mfem.FGMRESSolver())


def run_test():
    descriptor = json.loads(Path(registration.__file__).with_name(
        "ordinary_descriptor.json").read_text())
    classes = {entry["cxx_type"]: entry for entry in descriptor["classes"]}
    names = ("CGSolver", "GMRESSolver", "OperatorJacobiSmoother",
             "OperatorChebyshevSmoother", "SLISolver", "FGMRESSolver",
             "BiCGSTABSolver", "MINRESSolver", "NewtonSolver", "LBFGSSolver",
             "SLBQPOptimizer", "BlockILU", "DirectSubBlockSolver",
             "ProductSolver", "OrthoSolver")
    for name in names:
        entry = classes["mfem::" + name]
        assert any(method["cxx_name"] == "Mult" for method in entry["methods"]), name
    for name in ("IterativeSolver", "OptimizationSolver"):
        assert classes["mfem::" + name]["constructors"] == []
    for name in ("SLISolver", "CGSolver", "GMRESSolver", "FGMRESSolver",
                 "BiCGSTABSolver", "MINRESSolver", "NewtonSolver", "LBFGSSolver"):
        assert classes["mfem::" + name]["constructors"], name

    matrix = mfem.SparseMatrix(3)
    for i in range(3):
        matrix.Add(i, i, 1.0)
    matrix.Finalize()
    rhs = mfem.Vector(3)
    rhs.GetDataArray()[:] = [1.0, 2.0, 3.0]
    for name in ("SLISolver", "CGSolver", "GMRESSolver", "FGMRESSolver",
                 "BiCGSTABSolver", "MINRESSolver"):
        solver = make_fgmres() if name == "FGMRESSolver" else getattr(mfem, name)()
        solution = mfem.Vector(3)
        solution.Assign(0.0)
        assert solve(solver, matrix, rhs, solution), name
        np.testing.assert_allclose(solution.GetDataArray(), rhs.GetDataArray(),
                                   rtol=1e-11, atol=1e-12)
    assert solve.nopython_signatures and make_fgmres.nopython_signatures


if __name__ == "__main__":
    run_test()
