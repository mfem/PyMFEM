"""Compare the serial example-18 Numba director RHS with the Python version."""

from pathlib import Path
import sys

import numpy as np
import mfem.ser as mfem

registration = mfem.get_bridge_registration()

EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples"
BRIDGE_EXAMPLES_DIR = EXAMPLES_DIR / "bridge"
sys.path.insert(0, str(EXAMPLES_DIR))
sys.path.insert(0, str(BRIDGE_EXAMPLES_DIR))
import ex18 as bridge_ex18
from ex18_common import EulerInitialCondition, EulerMesh, DGHyperbolicConservationLaws


def make_integrator(dimension):
    flux = mfem.EulerFlux(dimension, 1.4)
    numerical_flux = mfem.RusanovFlux(flux)
    # Keep both proxy objects alive: HyperbolicFormIntegrator borrows the
    # numerical-flux instance, which in turn borrows EulerFlux.
    return flux, numerical_flux, mfem.HyperbolicFormIntegrator(numerical_flux, 1)


def run_test():
    mesh = EulerMesh("", 1)
    dimension = mesh.Dimension()
    equations = dimension + 2
    collection = mfem.DG_FECollection(1, dimension)
    space = mfem.FiniteElementSpace(mesh, collection, equations, mfem.Ordering.byNODES)
    solution = mfem.GridFunction(space)
    solution.ProjectCoefficient(EulerInitialCondition(1, 1.4, 1.0))

    reference_owners = make_integrator(dimension)
    reference = DGHyperbolicConservationLaws(
        space, reference_owners[2], preassembleWeakDivergence=True
    )
    bridge_owners = make_integrator(dimension)
    compiled = bridge_ex18.NumbaDGHyperbolicConservationLaws(
        space, bridge_owners[2], preassembleWeakDivergence=True
    )
    try:
        expected = mfem.Vector(solution.Size())
        actual = mfem.Vector(solution.Size())
        reference.Mult(solution, expected)
        compiled.Mult(solution, actual)
        assert np.max(np.abs(expected.GetDataArray() - actual.GetDataArray())) < 1.0e-11
        assert abs(reference.GetMaxCharSpeed() - compiled.GetMaxCharSpeed()) < 1.0e-11
    finally:
        compiled.close()


if __name__ == "__main__":
    run_test()
