"""Verify mfem::NonlinearForm calls from compiled Numba code."""
from numba import njit
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def zero_energy(form, values):
    form.Update()
    return form.GetEnergy(values), form.UseExternalIntegrators()


def run_test():
    mesh = mfem.Mesh(1, 1, "QUADRILATERAL")
    fec = mfem.H1_FECollection(1, mesh.Dimension())
    space = mfem.FiniteElementSpace(mesh, fec)
    form = mfem.NonlinearForm(space)
    values = mfem.Vector(space.GetVSize())
    values.Assign(0.0)
    assert zero_energy(form, values) == (0.0, None)


if __name__ == "__main__":
    run_test()
