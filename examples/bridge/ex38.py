"""Serial example 38 using direct NSB LinearFormIntegrator directors.

The existing ``examples/ex38.py`` supplies the moment-fitting rules and the
``mfem.jit.scalar`` level-set and integrand coefficients. This version moves
the two element quadrature loops into compiled direct MFEM directors.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from numba import types
import numba_swig_bridge as nsb
import mfem.ser as mfem


registration = mfem.get_bridge_registration()
EXAMPLES = Path(__file__).resolve().parents[1]


def _load_reference():
    spec = importlib.util.spec_from_file_location("_pymfem_ex38_reference",
                                                  EXAMPLES / "ex38.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.visualization = False
    return module


reference = _load_reference()
ARRAY2D = types.float64[::1, :]
ARRAY3D = types.Array(types.float64, 3, "A")
STATE = {
    "coefficient": mfem.Coefficient,
    "weights": ARRAY2D,
    "coordinates": ARRAY3D,
    "shape": mfem.Vector,
    "point": mfem.IntegrationPoint,
    "npoints": types.intc,
}


def _initialize(integrator, coefficient, weights, coordinates):
    integrator.coefficient = coefficient
    integrator.weights = np.asfortranarray(weights, dtype=np.float64)
    integrator.coordinates = np.asarray(coordinates, dtype=np.float64)
    integrator.shape = mfem.Vector()
    integrator.point = mfem.IntegrationPoint()
    integrator.npoints = weights.shape[0]


@nsb.director(state=STATE, fallback="silent")
class SurfaceLFIntegrator(mfem.LinearFormIntegrator):
    def __init__(self, coefficient, weights, coordinates):
        _initialize(self, coefficient, weights, coordinates)
        super().__init__()

    @nsb.override(
        stateaccess="coefficient, weights, coordinates, shape, point, npoints"
    )
    def AssembleRHSElementVect(self, element, transformation, element_vector):
        dof = element.GetDof()
        self.shape.SetSize(dof)
        element_vector.SetSize(dof)
        output = element_vector.GetDataArray()
        shape = self.shape.GetDataArray()
        for j in range(dof):
            output[j] = 0.0
        element_number = transformation.ElementNo
        for i in range(self.npoints):
            self.point.x = self.coordinates[i, element_number, 0]
            self.point.y = self.coordinates[i, element_number, 1]
            self.point.z = self.coordinates[i, element_number, 2]
            transformation.SetIntPoint(self.point)
            element.CalcShape(self.point, self.shape)
            scale = (self.weights[i, element_number] * transformation.Weight()
                     * self.coefficient.Eval(transformation, self.point))
            for j in range(dof):
                output[j] += scale * shape[j]


@nsb.director(state=STATE, fallback="silent")
class SubdomainLFIntegrator(mfem.LinearFormIntegrator):
    def __init__(self, coefficient, weights, coordinates):
        _initialize(self, coefficient, weights, coordinates)
        super().__init__()

    @nsb.override(
        stateaccess="coefficient, weights, coordinates, shape, point, npoints"
    )
    def AssembleRHSElementVect(self, element, transformation, element_vector):
        dof = element.GetDof()
        self.shape.SetSize(dof)
        element_vector.SetSize(dof)
        output = element_vector.GetDataArray()
        shape = self.shape.GetDataArray()
        for j in range(dof):
            output[j] = 0.0
        element_number = transformation.ElementNo
        for i in range(self.npoints):
            self.point.x = self.coordinates[i, element_number, 0]
            self.point.y = self.coordinates[i, element_number, 1]
            self.point.z = self.coordinates[i, element_number, 2]
            transformation.SetIntPoint(self.point)
            element.CalcPhysShape(transformation, self.shape)
            scale = (self.weights[i, element_number] * transformation.Weight()
                     * self.coefficient.Eval(transformation, self.point))
            for j in range(dof):
                output[j] += scale * shape[j]


def _integration_type(name):
    normalized = name.lower()
    choices = {
        "volumetric1d": reference.IntegrationType.Volumetric1D,
        "surface2d": reference.IntegrationType.Surface2D,
        "volumetric2d": reference.IntegrationType.Volumetric2D,
        "surface3d": reference.IntegrationType.Surface3D,
        "volumetric3d": reference.IntegrationType.Volumetric3D,
    }
    return choices[normalized]


def _make_mesh(integration_type):
    if integration_type == reference.IntegrationType.Volumetric1D:
        return mfem.Mesh(str(EXAMPLES.parent / "data" / "inline-segment.mesh"))
    if integration_type in (reference.IntegrationType.Surface2D,
                             reference.IntegrationType.Volumetric2D):
        mesh = mfem.Mesh(2, 4, 1, 0, 2)
        mesh.AddVertex(-1.6, -1.6)
        mesh.AddVertex(1.6, -1.6)
        mesh.AddVertex(1.6, 1.6)
        mesh.AddVertex(-1.6, 1.6)
        mesh.AddQuad(0, 1, 2, 3)
        mesh.FinalizeQuadMesh(1, 0, True)
        return mesh
    mesh = mfem.Mesh(3, 8, 1, 0, 3)
    mesh.AddVertex(-1.6, -1.6, -1.6)
    mesh.AddVertex(1.6, -1.6, -1.6)
    mesh.AddVertex(1.6, 1.6, -1.6)
    mesh.AddVertex(-1.6, 1.6, -1.6)
    mesh.AddVertex(-1.6, -1.6, 1.6)
    mesh.AddVertex(1.6, -1.6, 1.6)
    mesh.AddVertex(1.6, 1.6, 1.6)
    mesh.AddVertex(-1.6, 1.6, 1.6)
    mesh.AddHex(0, 1, 2, 3, 4, 5, 6, 7)
    mesh.FinalizeHexMesh(1, 0, True)
    return mesh


def _quadrature_state(rule, elements, surface):
    count = rule.GetNPoints()
    coordinates = np.empty((count, elements, 3), dtype=np.float64)
    weights = np.empty((count, elements), dtype=np.float64, order="F")
    raw_weights = rule.Weights.GetDataArray()
    surface_weights = rule.SurfaceWeights.GetDataArray() if surface else None
    for element in range(elements):
        for i in range(count):
            point = rule.IntPoint(i)
            coordinates[i, element, 0] = point.x
            coordinates[i, element, 1] = point.y
            coordinates[i, element, 2] = point.z
            if rule.dim == 1:
                offset = 0 if surface else 2 * i
                coordinates[i, element, 0] = raw_weights[offset, element]
                weights[i, element] = raw_weights[offset + 1, element]
            else:
                weights[i, element] = raw_weights[i, element]
                if surface:
                    weights[i, element] *= surface_weights[i, element]
    return weights, coordinates


def run(backend="bridge", integration_type="surface2d", ref_levels=3, order=2):
    selected = _integration_type(integration_type)
    reference.itype = selected
    mesh = _make_mesh(selected)
    for _ in range(ref_levels):
        mesh.UniformRefinement()
    collection = mfem.H1_FECollection(1, mesh.Dimension())
    space = mfem.FiniteElementSpace(mesh, collection)
    levelset = reference.make_lvlset_coeff()
    coefficient = reference.make_integrand_coeff()
    surface_rule = reference.SIntegrationRule(order, levelset, 2, mesh)
    surface_state = _quadrature_state(surface_rule, mesh.GetNE(), True)
    volumetric = selected in (
        reference.IntegrationType.Volumetric1D,
        reference.IntegrationType.Volumetric2D,
        reference.IntegrationType.Volumetric3D,
    )
    volume_rule = None
    volume_state = None
    if volumetric:
        volume_rule = reference.CIntegrationRule(order, levelset, 2, mesh)
        volume_state = _quadrature_state(volume_rule, mesh.GetNE(), False)

    surface = mfem.LinearForm(space)
    if backend == "bridge":
        surface_integrator = SurfaceLFIntegrator(
            coefficient, surface_state[0], surface_state[1])
    else:
        surface_integrator = reference.SurfaceLFIntegrator(
            coefficient, levelset, surface_rule)
    surface.AddDomainIntegrator(surface_integrator)
    if backend == "bridge":
        assert surface_integrator.thisown is False
        surface_integrator.close()
    start = perf_counter()
    surface.Assemble()
    surface_seconds = perf_counter() - start

    volume = None
    volume_seconds = None
    if volumetric:
        volume = mfem.LinearForm(space)
        if backend == "bridge":
            volume_integrator = SubdomainLFIntegrator(
                coefficient, volume_state[0], volume_state[1])
        else:
            volume_integrator = reference.SubdomainLFIntegrator(
                coefficient, levelset, volume_rule)
        volume.AddDomainIntegrator(volume_integrator)
        if backend == "bridge":
            assert volume_integrator.thisown is False
            volume_integrator.close()
        start = perf_counter()
        volume.Assemble()
        volume_seconds = perf_counter() - start

    return {
        "surface": surface.GetDataArray().copy(),
        "surface_sum": surface.Sum(),
        "surface_assembly": surface_seconds,
        "volume": None if volume is None else volume.GetDataArray().copy(),
        "volume_sum": None if volume is None else volume.Sum(),
        "volume_assembly": volume_seconds,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--integrationtype", default="surface2d",
                        choices=("volumetric1d", "surface2d", "volumetric2d",
                                 "surface3d", "volumetric3d"))
    parser.add_argument("-r", "--refine", default=3, type=int)
    parser.add_argument("-o", "--order", default=2, type=int)
    parser.add_argument("--backend", choices=("bridge", "python"), default="bridge")
    parser.add_argument("--compare", action="store_true")
    args = parser.parse_args()
    result = run(args.backend, args.integrationtype, args.refine, args.order)
    print(json.dumps({key: value for key, value in result.items()
                      if not isinstance(value, np.ndarray)}))
    if args.compare:
        expected = run("python", args.integrationtype, args.refine, args.order)
        np.testing.assert_allclose(result["surface"], expected["surface"],
                                   rtol=1e-11, atol=1e-12)
        if result["volume"] is not None:
            np.testing.assert_allclose(result["volume"], expected["volume"],
                                       rtol=1e-11, atol=1e-12)
        print("Python and NSB linear forms agree.")


if __name__ == "__main__":
    if not hasattr(mfem, "MomentFittingIntRules"):
        raise RuntimeError("MFEM must be built with LAPACK for example 38")
    main()
