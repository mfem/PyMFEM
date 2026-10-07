"""Run MFEM example 18 with a Numba-compiled volume callback.

This serial demonstration retains the ordinary PyMFEM setup and ODE solver.
Only the per-element volume contribution in ``Mult`` is a native Numba
Director callback. It supports the preassembled weak-divergence mode and
uniform-order finite-element spaces, which are the default example-18 path.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import numba_swig_bridge as nsb
from numba import njit, types
import mfem.ser as mfem

# Install the generated serial ordinary-method and director registrations.
registration = mfem.get_bridge_registration()

EXAMPLES_DIR = Path(__file__).resolve().parents[1]
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))
import ex18 as reference_ex18


@nsb.director(
    state={
        "vfes": mfem.FiniteElementSpace,
        "form_integrator": mfem.HyperbolicFormIntegrator,
        "nonlinear_form": mfem.NonlinearForm,
        "z": mfem.Vector,
        "invmass": mfem.DenseTensor,
        "weakdiv": mfem.DenseTensor,
        "vdofs": mfem.intArray,
        "num_equations": types.intc,
        "dimension": types.intc,
        "element_count": types.intc,
        "element_dofs": types.intc,
        "x_data": types.float64[:, ::1],
        "z_data": types.float64[:, ::1],
        "y_data": types.float64[:, ::1],
        "flux_data": types.float64[:, ::1],
        "current_flux_data": types.float64[:, ::1],
        "state_data": types.float64[::1],
    },
    # Automatic director selection also exposes non-pure MFEM virtuals. This
    # example implements only Mult; the rest intentionally use MFEM defaults.
    fallback="silent",
)
class NumbaDGHyperbolicConservationLaws(mfem.TimeDependentOperator):
    """Example-18 operator with NumPy work arrays inside a Numba director.

    Element matrices stay in C++-owned tensors. NumPy work buffers live in
    StateBlock and are reused across calls; frame-owned MFEM views expose the
    state and flux buffers for ``ComputeFlux``. This serial operator assumes
    calls on each operator are neither concurrent nor reentrant.
    """

    def __init__(self, vfes, form_integrator, preassembleWeakDivergence=True):
        if not preassembleWeakDivergence:
            raise ValueError("bridge ex18 currently requires preassembled weak divergence")
        super().__init__(vfes.GetTrueVSize())
        element_count = vfes.GetNE()
        element_dofs = vfes.GetFE(0).GetDof()
        for element in range(1, element_count):
            if vfes.GetFE(element).GetDof() != element_dofs:
                raise ValueError("bridge ex18 requires a uniform finite-element order")

        self.vfes = vfes
        self.form_integrator = form_integrator
        self.num_equations = form_integrator.num_equations
        self.dimension = vfes.GetMesh().SpaceDimension()
        self.element_count = element_count
        self.element_dofs = element_dofs
        # Fixed-size reusable buffers for this uniform-order serial operator.
        shape = (self.num_equations, element_dofs)
        self.x_data = np.empty(shape, dtype=np.float64)
        self.z_data = np.empty(shape, dtype=np.float64)
        self.y_data = np.empty(shape, dtype=np.float64)
        self.flux_data = np.empty(
            (self.num_equations, self.dimension * element_dofs), dtype=np.float64
        )
        self.state_data = np.empty(self.num_equations, dtype=np.float64)
        self.current_flux_data = np.empty(
            (self.dimension, self.num_equations), dtype=np.float64
        )
        self.z = mfem.Vector(vfes.GetTrueVSize())
        self.vdofs = mfem.intArray()
        self.invmass = mfem.DenseTensor(element_dofs, element_dofs, element_count)
        self.weakdiv = mfem.DenseTensor(
            element_dofs, element_dofs * self.dimension, element_count
        )

        self._assemble_element_data()
        self._validate_vdofs()
        nonlinear_form = mfem.NonlinearForm(vfes)
        nonlinear_form.AddInteriorFaceIntegrator(form_integrator)
        nonlinear_form.UseExternalIntegrators()
        self.nonlinear_form = nonlinear_form

    def _validate_vdofs(self):
        """Require unsigned DG vdofs before using direct data-view indexing."""
        for element in range(self.element_count):
            self.vfes.GetElementVDofs(element, self.vdofs)
            if any(index < 0 for index in self.vdofs.ToList()):
                raise ValueError(
                    "bridge ex18 data-view gather/scatter requires unsigned element vdofs"
                )

    def _assemble_element_data(self):
        """Preassemble Python-side matrices into C++ tensors for the callback."""
        inv_mass = mfem.InverseIntegrator(mfem.MassIntegrator())
        weak_div = mfem.TransposeIntegrator(mfem.GradientIntegrator())
        mass = mfem.DenseMatrix()
        weak_by_nodes = mfem.DenseMatrix()
        invmass_data = self.invmass.GetDataArray()
        weakdiv_data = self.weakdiv.GetDataArray()
        for element in range(self.element_count):
            finite_element = self.vfes.GetFE(element)
            transformation = self.vfes.GetElementTransformation(element)
            inv_mass.AssembleElementMatrix(finite_element, transformation, mass)
            invmass_data[element, :, :] = mass.GetDataArray()

            weak_by_nodes.SetSize(
                self.element_dofs, self.element_dofs * self.dimension
            )
            weak_div.AssembleElementMatrix2(
                finite_element, finite_element, transformation, weak_by_nodes
            )
            source = weak_by_nodes.GetDataArray()
            # Preserve ex18_common.py's ByDim-to-ByNodes column ordering.
            for dof in range(self.element_dofs):
                for direction in range(self.dimension):
                    destination_column = dof * self.dimension + direction
                    source_column = direction * self.element_dofs + dof
                    weakdiv_data[element, :, destination_column] = source[:, source_column]

    def GetMaxCharSpeed(self):
        return self.form_integrator.GetMaxCharSpeed()

    @nsb.override(stateaccess="*")
    def Mult(self, x, y):
        self.form_integrator.ResetMaxCharSpeed()
        self.nonlinear_form.Mult(x, self.z)

        inverse_mass = self.invmass.GetDataArray()
        weak_divergence = self.weakdiv.GetDataArray()
        input_data = x.GetDataArray()
        auxiliary_data = self.z.GetDataArray()
        output_data = y.GetDataArray()
        local_size = self.element_dofs * self.num_equations

        # We make local variables. It looks like this matters for speed (order of 5-7%)
        x_data = self.x_data
        z_data = self.z_data
        y_data = self.y_data
        flux_data = self.flux_data
        state_data = self.state_data
        current_flux_data = self.current_flux_data
        x_data_flat = x_data.reshape(local_size)
        z_data_flat = z_data.reshape(local_size)
        y_data_flat = y_data.reshape(local_size)
        current_state = mfem.Vector(
            nsb.as_cpointer(state_data), self.num_equations
        )
        current_flux = mfem.DenseMatrix(
            nsb.as_cpointer(current_flux_data),
            self.num_equations,
            self.dimension,
        )
        flux_function = self.form_integrator.GetFluxFunction()
        # The generic bridge maps C++ pointer returns to Optional[T]. MFEM
        # creates this helper with the integrator, but retain the native null
        # contract rather than dereferencing a nullable borrowed object.
        flux_function = nsb.require_not_none(flux_function)

        for element in range(self.element_count):
            transformation = nsb.require_not_none(
                self.vfes.GetElementTransformation(element)
            )
            # The nullable borrowed DofTransformation is irrelevant here; the
            # mutable Array<int>& receives the element vector dofs.
            self.vfes.GetElementVDofs(element, self.vdofs)
            vdof_data = self.vdofs.GetDataArray()
            # Direct gather avoids fancy-index temporaries and slice broadcasting.
            for local in range(local_size):
                x_data_flat[local] = input_data[vdof_data[local]]
                z_data_flat[local] = auxiliary_data[vdof_data[local]]

            for node in range(self.element_dofs):
                for equation in range(self.num_equations):
                    state_data[equation] = x_data[equation, node]
                flux_function.ComputeFlux(
                    current_state, transformation, current_flux
                )
                start = self.dimension * node
                for equation in range(self.num_equations):
                    for direction in range(self.dimension):
                        flux_data[equation, start + direction] = (
                            current_flux_data[direction, equation]
                        )

            # This is the transposed NumPy form of the original MFEM calls:
            # AddMult_a_ABt(1.0, weakdiv, flux, current_zmat) and
            # Mult(invmass, current_zmat, current_ymat).
            # A DenseTensor is a single contiguous MFEM buffer whose 2-D
            # element slices use MFEM's column-major matrix layout.  The
            # complete 3-D view is therefore strided, but each slice is
            # Fortran-contiguous.  These calls only assert that slice layout;
            # they do not copy its data.
            weak_divergence_element = np.asfortranarray(
                weak_divergence[element]
            )
            inverse_mass_element = np.asfortranarray(inverse_mass[element])
            z_data += flux_data @ weak_divergence_element.T
            y_data[:] = z_data @ inverse_mass_element.T
            for local in range(local_size):
                output_data[vdof_data[local]] = y_data_flat[local]




@njit
def _apply_numba_operator(operator, source, destination):
    operator.Mult(source, destination)


def _python_mult(self, source, destination):
    """Route direct Python calls through C++ virtual dispatch as well."""
    _apply_numba_operator(self, source, destination)


# The director decorator has already captured the compiled ``Mult`` callback.
# Replacing the Python-facing method lets ex18's CFL preflight use that same
# native callback instead of attempting to execute the Numba-only body.
NumbaDGHyperbolicConservationLaws.Mult = _python_mult


def run(**kwargs):
    """Run the standard serial example using the Numba director operator."""
    original = reference_ex18.DGHyperbolicConservationLaws
    reference_ex18.DGHyperbolicConservationLaws = NumbaDGHyperbolicConservationLaws
    try:
        return reference_ex18.run(preassembleWeakDiv=True, **kwargs)
    finally:
        reference_ex18.DGHyperbolicConservationLaws = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-p", "--problem", type=int, default=1)
    parser.add_argument("-r", "--refine", type=int, default=0)
    parser.add_argument("-o", "--order", type=int, default=1)
    parser.add_argument("-s", "--ode-solver", dest="ode_solver", type=int, default=4)
    parser.add_argument("-tf", "--t-final", dest="t_final", type=float, default=0.02)
    parser.add_argument("-c", "--cfl-number", dest="cfl", type=float, default=0.3)
    parser.add_argument("-vs", "--visualization-steps", dest="vis_steps", type=int, default=50)
    parser.add_argument("-m", "--mesh", default="")
    parser.add_argument("--visualization", action="store_true")
    args = parser.parse_args()
    run(
        problem=args.problem,
        ref_levels=args.refine,
        order=args.order,
        ode_solver_type=args.ode_solver,
        t_final=args.t_final,
        cfl=args.cfl,
        visualization=args.visualization,
        vis_steps=args.vis_steps,
        meshfile=args.mesh,
    )


if __name__ == "__main__":
    main()
