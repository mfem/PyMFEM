"""GridFunction inheritance, methods, and state mutation; use -p for MPI."""

import sys

import numpy as np
from numba import njit, typeof
import nsb_rt as nsb

if "-p" in sys.argv:
    from mpi4py import MPI
    import mfem.par as mfem
else:
    import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def mutate(field, other, output):
    field[0] = 3.0
    field += other
    field -= 1.0
    field *= 2.0
    field /= 2.0
    mfem.Add(field, other, output)
    space = nsb.require_not_none(field.FESpace())
    return field.Size(), space.GetVSize(), field * other


@njit
def true_dofs(field, values):
    field.GetTrueDofs(values)
    field.SetFromTrueDofs(values)
    return values.Size()


@njit
def parallel_space(field):
    return nsb.require_not_none(field.ParFESpace()).GetVSize()


@njit
def invoke(operator, source, destination):
    operator.Mult(source, destination)


@njit
def construct(space):
    field = mfem.GridFunction(space)
    for i in range(field.Size()):
        field[i] = 1.0
    field += 2.0
    return nsb.own_this(field)


@njit
def construct_parallel(space):
    field = mfem.ParGridFunction(space)
    for i in range(field.Size()):
        field[i] = 1.0
    field += 2.0
    return nsb.own_this(field)


@njit
def allocated_true_dofs(field):
    return nsb.own_this(nsb.require_not_none(field.GetTrueDofs()))


def exercise(field_type, space):
    field, other, output = (field_type(space) for _ in range(3))
    field.Assign(2.0)
    other.Assign(2.0)
    # Explicitly prepare native host indexing; the bridge adds no synchronization.
    field.HostReadWrite()
    other.HostRead()
    size = space.GetVSize()
    assert mutate(field, other, output) == (size, size, 2.0 * (3.0 * size + 1.0))
    np.testing.assert_allclose(field.GetDataArray(), [4.0] + [3.0] * (size - 1))
    np.testing.assert_allclose(output.GetDataArray(), [6.0] + [5.0] * (size - 1))
    values = mfem.Vector()
    assert true_dofs(field, values) == space.GetTrueVSize()
    before = field.GetDataArray().copy()

    @nsb.director(state={"field": field_type}, fallback="silent")
    class Mutator(mfem.Operator):
        def __init__(self, value):
            self.field = value
            super().__init__(value.Size())

        @nsb.override
        def Mult(self, source, destination):
            self.field += 2.0
            self.field -= 1.0
            self.field *= 2.0
            self.field /= 2.0
            self.field.GetTrueDofs(destination)

    mutator = Mutator(field)
    try:
        owned_before = field.thisown
        invoke(mutator, other, values)
        assert field.thisown == owned_before
        assert typeof(mutator.field).cxx_type == typeof(field).cxx_type
        np.testing.assert_allclose(field.GetDataArray(), before + 1.0)
    finally:
        mutator.close()
    if field_type is getattr(mfem, "ParGridFunction", None):
        assert parallel_space(field) == size
        allocated = allocated_true_dofs(field)
        assert allocated.thisown
        np.testing.assert_allclose(allocated.GetDataArray(), values.GetDataArray())


def run_test():
    mesh = mfem.Mesh(2, 1, "QUADRILATERAL")
    fec = mfem.H1_FECollection(1, mesh.Dimension())
    space = mfem.FiniteElementSpace(mesh, fec)
    exercise(mfem.GridFunction, space)
    owned = construct(space)
    assert owned.thisown
    np.testing.assert_allclose(owned.GetDataArray(), 3.0)
    if "-p" in sys.argv:
        pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
        pspace = mfem.ParFiniteElementSpace(pmesh, fec)
        exercise(mfem.ParGridFunction, pspace)
        powned = construct_parallel(pspace)
        assert powned.thisown
        np.testing.assert_allclose(powned.GetDataArray(), 3.0)
    print("GridFunction bridges: Normal End")


if __name__ == "__main__":
    run_test()
