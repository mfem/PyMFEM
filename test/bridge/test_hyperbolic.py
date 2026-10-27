"""Verify FluxFunction and HyperbolicFormIntegrator ordinary bridges."""
from numba import njit
import mfem.ser as mfem

registration = mfem.get_bridge_registration()


@njit
def burgers_flux(integrator, state, transformation, flux):
    integrator.ResetMaxCharSpeed()
    speed = integrator.GetFluxFunction().ComputeFlux(state, transformation, flux)
    return speed, integrator.GetMaxCharSpeed()


def run_test():
    mesh = mfem.Mesh(1, 1, "QUADRILATERAL")
    transformation = mesh.GetElementTransformation(0)
    state = mfem.Vector(1)
    state[0] = 3.0
    flux = mfem.DenseMatrix(1, 2)
    # HyperbolicFormIntegrator borrows its NumericalFlux, which in turn
    # borrows its FluxFunction. Keep both Python owners alive for the test.
    flux_function = mfem.BurgersFlux(2)
    numerical_flux = mfem.RusanovFlux(flux_function)
    integrator = mfem.HyperbolicFormIntegrator(numerical_flux, 1)

    assert burgers_flux(integrator, state, transformation, flux) == (3.0, 0.0)
    assert flux[0, 0] == 4.5
    assert flux[0, 1] == 4.5


if __name__ == "__main__":
    run_test()
