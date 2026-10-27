"""Verify finite-element-space calls and callback-scoped borrowed returns."""
from numba import njit
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def element_metadata(space):
    element = space.GetFE(0)
    transformation = space.GetElementTransformation(0)
    return (space.GetVSize(), element.GetDof(), transformation.GetDimension(),
            transformation.GetSpaceDim())


@njit
def element_vdofs(space, vdofs):
    # GetElementVDofs returns an optional borrowed DofTransformation.  The
    # result can be null; the Array<int>& is still populated in either case.
    transformation = space.GetElementVDofs(0, vdofs)
    if transformation is None:
        return vdofs.Size()
    return vdofs.Size() + transformation.Size()


def run_test():
    mesh = mfem.Mesh(1, 1, "QUADRILATERAL")
    fec = mfem.H1_FECollection(1, mesh.Dimension())
    space = mfem.FiniteElementSpace(mesh, fec)
    assert element_metadata(space) == (4, 4, 2, 2)
    vdofs = mfem.intArray()
    assert element_vdofs(space, vdofs) == 4
    assert sorted(vdofs.ToList()) == [0, 1, 2, 3]


if __name__ == "__main__":
    run_test()
