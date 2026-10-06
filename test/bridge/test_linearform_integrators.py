"""Direct LinearFormIntegrator directors through real MFEM assembly."""

import gc

import numpy as np
from numba import types
import numba_swig_bridge as nsb
import mfem.ser as mfem


registration = mfem.get_bridge_registration()
ARRAY2D = types.float64[::1, :]


def _initialize(integrator, coefficient, rule, weights):
    integrator.coefficient = coefficient
    integrator.rule = rule
    integrator.weights = np.asfortranarray(weights, dtype=np.float64)
    integrator.shape = mfem.Vector()
    integrator.point = mfem.IntegrationPoint()


STATE = {
    "coefficient": mfem.Coefficient,
    "rule": mfem.IntegrationRule,
    "weights": ARRAY2D,
    "shape": mfem.Vector,
    "point": mfem.IntegrationPoint,
}


@nsb.director(state=STATE, fallback="silent")
class SurfaceIntegrator(mfem.LinearFormIntegrator):
    def __init__(self, coefficient, rule, weights):
        _initialize(self, coefficient, rule, weights)
        super().__init__()

    @nsb.override(stateaccess="coefficient, rule, weights, shape, point")
    def AssembleRHSElementVect(self, element, transformation, element_vector):
        dof = element.GetDof()
        self.shape.SetSize(dof)
        element_vector.SetSize(dof)
        output = element_vector.GetDataArray()
        for j in range(dof):
            output[j] = 0.0
        shape = self.shape.GetDataArray()
        element_number = transformation.ElementNo
        for i in range(self.rule.GetNPoints()):
            source = self.rule.IntPoint(i)
            self.point.x = source.x
            self.point.y = source.y
            self.point.z = source.z
            transformation.SetIntPoint(self.point)
            element.CalcShape(self.point, self.shape)
            scale = (self.weights[i, element_number] * transformation.Weight()
                     * self.coefficient.Eval(transformation, self.point))
            for j in range(dof):
                output[j] += scale * shape[j]


@nsb.director(state=STATE, fallback="silent")
class SubdomainIntegrator(mfem.LinearFormIntegrator):
    def __init__(self, coefficient, rule, weights):
        _initialize(self, coefficient, rule, weights)
        super().__init__()

    @nsb.override(stateaccess="coefficient, rule, weights, shape, point")
    def AssembleRHSElementVect(self, element, transformation, element_vector):
        dof = element.GetDof()
        self.shape.SetSize(dof)
        element_vector.SetSize(dof)
        output = element_vector.GetDataArray()
        for j in range(dof):
            output[j] = 0.0
        shape = self.shape.GetDataArray()
        element_number = transformation.ElementNo
        for i in range(self.rule.GetNPoints()):
            source = self.rule.IntPoint(i)
            self.point.x = source.x
            self.point.y = source.y
            self.point.z = source.z
            transformation.SetIntPoint(self.point)
            element.CalcPhysShape(transformation, self.shape)
            scale = (self.weights[i, element_number] * transformation.Weight()
                     * self.coefficient.Eval(transformation, self.point))
            for j in range(dof):
                output[j] += scale * shape[j]


class PythonIntegrator(mfem.PyLinearFormIntegrator):
    def __init__(self, coefficient, rule, weights, physical):
        super().__init__()
        self.coefficient = coefficient
        self.rule = rule
        self.weights = weights
        self.physical = physical
        self.shape = mfem.Vector()

    def AssembleRHSElementVect(self, element, transformation, element_vector):
        dof = element.GetDof()
        self.shape.SetSize(dof)
        element_vector.SetSize(dof)
        element_vector.Assign(0.0)
        element_number = transformation.ElementNo
        for i in range(self.rule.GetNPoints()):
            point = self.rule.IntPoint(i)
            transformation.SetIntPoint(point)
            if self.physical:
                element.CalcPhysShape(transformation, self.shape)
            else:
                element.CalcShape(point, self.shape)
            scale = (self.weights[i, element_number] * transformation.Weight()
                     * self.coefficient.Eval(transformation, point))
            mfem.add_vector(element_vector, scale, self.shape, element_vector)


def _assemble(space, integrator):
    form = mfem.LinearForm(space)
    form.AddDomainIntegrator(integrator)
    assert integrator.thisown is False
    state = integrator._nsb_director_state if hasattr(integrator, "_nsb_director_state") else None
    if state is not None:
        integrator.close()
    form.Assemble()
    values = form.GetDataArray().copy()
    del form
    gc.collect()
    if state is not None:
        assert not state.IsAlive()
    return values


def run_test():
    mesh = mfem.Mesh(2, 2, "QUADRILATERAL")
    mesh.UniformRefinement()
    fec = mfem.H1_FECollection(2, 2)
    space = mfem.FiniteElementSpace(mesh, fec)
    coefficient = mfem.ConstantCoefficient(2.5)
    rules = mfem.IntegrationRules(0, mfem.Quadrature1D.GaussLegendre)
    rule = rules.Get(mfem.Geometry.SQUARE, 4)
    weights = np.empty((rule.GetNPoints(), mesh.GetNE()), order="F")
    for element in range(mesh.GetNE()):
        for i in range(rule.GetNPoints()):
            weights[i, element] = rule.IntPoint(i).weight

    surface = _assemble(space, SurfaceIntegrator(coefficient, rule, weights))
    surface_python = _assemble(
        space, PythonIntegrator(coefficient, rule, weights, physical=False))
    subdomain = _assemble(space, SubdomainIntegrator(coefficient, rule, weights))
    subdomain_python = _assemble(
        space, PythonIntegrator(coefficient, rule, weights, physical=True))

    native_form = mfem.LinearForm(space)
    native_form.AddDomainIntegrator(mfem.DomainLFIntegrator(coefficient, rule))
    native_form.Assemble()
    native = native_form.GetDataArray().copy()

    np.testing.assert_allclose(surface, surface_python, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(subdomain, subdomain_python, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(subdomain, native, rtol=1e-13, atol=1e-13)


if __name__ == "__main__":
    run_test()
