"""Rank-local ordinary/director/MPI bridge acceptance; run on two MPI ranks."""

import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from numba import njit, types
import nsb_rt as nsb
import mfem.par as mfem

registration = mfem.get_bridge_registration()


@njit
def owned_vector(size):
    return nsb.own_this(mfem.Vector(size))


@njit
def comm_info(comm):
    return comm.rank, comm.size


@njit
def make_hypre(comm):
    solver = mfem.HyprePCG(comm)
    solver.SetTol(1e-12)
    solver.SetMaxIter(8)
    solver.SetPrintLevel(0)
    return nsb.own_this(solver)


@njit
def inspect_comm(matrix):
    comm = matrix.GetComm()
    solver = mfem.CGSolver(comm)
    solver.SetPrintLevel(0)
    return comm.rank, comm.size


@njit
def apply(operator, source, destination):
    operator.Mult(source, destination)


@nsb.director(state={"fail": types.boolean}, fallback="silent")
class Inner(mfem.Operator):
    def __init__(self):
        super().__init__(3)
        self.fail = False

    @nsb.override
    def Mult(self, source, destination):
        if self.fail:
            raise ValueError("rank-local callback failure")
        x = source.GetDataArray()
        y = destination.GetDataArray()
        for i in range(x.size):
            y[i] = 2.0 * x[i]


@nsb.director(state={"inner": mfem.Operator}, fallback="silent")
class Outer(mfem.Operator):
    def __init__(self, inner):
        super().__init__(3)
        self.inner = inner

    @nsb.override
    def Mult(self, source, destination):
        self.inner.Mult(source, destination)


def run_test():
    comm = MPI.COMM_WORLD
    assert comm.size >= 2
    assert comm_info(comm) == (comm.rank, comm.size)
    hypre_solver = make_hypre(comm)
    assert isinstance(hypre_solver, mfem.HyprePCG)
    del hypre_solver
    descriptor = json.loads(Path(registration.__file__).with_name(
        "ordinary_descriptor.json").read_text())
    classes = {entry["cxx_type"]: entry for entry in descriptor["classes"]}
    for name in ("HypreSmoother", "HypreSolver", "HypreTriSolve", "HyprePCG",
                 "HypreGMRES", "HypreFGMRES", "HypreIdentity", "HypreDiagScale",
                 "HypreParaSails", "HypreEuclid", "HypreILU", "HypreBoomerAMG",
                 "HypreAMS", "HypreADS", "HypreLOBPCG", "HypreAME"):
        assert "mfem::" + name in classes, name
    for name in ("HypreLOBPCG", "HypreAME"):
        assert any("StealEigenvectors" in item["signature"]
                   for item in classes["mfem::" + name]["skipped_methods"])

    source = owned_vector(3)
    destination = mfem.Vector(3)
    source.Assign(comm.rank + 1.0)
    inner = Inner()
    outer = Outer(inner)
    try:
        apply(outer, source, destination)
        np.testing.assert_array_equal(destination.GetDataArray(), 2 * source.GetDataArray())
        inner.fail = True
        try:
            apply(inner, source, destination)
        except RuntimeError as error:
            # Ordinary adapters report callback failure through their status ABI.
            assert "ordinary-object bridge adapter failed" in str(error), str(error)
        else:
            raise AssertionError("callback error was not propagated")
        inner.fail = False
        apply(outer, source, destination)
    finally:
        outer.close()
        inner.close()

    mesh = mfem.Mesh.MakeCartesian2D(2, 2, mfem.Element.QUADRILATERAL)
    pmesh = mfem.ParMesh(comm, mesh)
    fec = mfem.H1_FECollection(1, 2)
    space = mfem.ParFiniteElementSpace(pmesh, fec)
    form = mfem.ParBilinearForm(space)
    form.AddDomainIntegrator(mfem.MassIntegrator())
    form.Assemble()
    form.Finalize()
    matrix = form.ParallelAssemble()
    assert inspect_comm(matrix) == (comm.rank, comm.size)
    del matrix, form, space, fec, pmesh, mesh, outer, inner, destination, source
    comm.Barrier()
    if comm.rank == 0:
        print("Parallel foundation passed on", comm.size, "ranks", flush=True)


if __name__ == "__main__":
    run_test()
