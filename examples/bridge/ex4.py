"""Serial example 4 using direct NSB VectorCoefficient directors.

Run --compare to check Python callbacks and repeat the solve with warmed
callbacks. Assembly and projection timings exclude imports and construction.
"""
from time import perf_counter
IMPORT_STARTED = perf_counter()

import argparse
import json
from pathlib import Path
import numpy as np
from numpy import sin, cos
from numba import types
import numba_swig_bridge as nsb
import mfem.ser as mfem
from mfem.ser import intArray

registration = mfem.get_bridge_registration()
kappa = np.pi
CALLBACK_COMPILE_STARTED = perf_counter()


@nsb.director(state={"point": mfem.Vector, "dimension": types.intc,
                     "scale": types.float64}, fallback="silent")
class NumbaCoefficient(mfem.VectorCoefficient):
    def __init__(self, dimension, scale):
        super().__init__(dimension)
        self.point = mfem.Vector(dimension)
        self.dimension = dimension
        self.scale = scale

    @nsb.override
    def Eval(self, output, transformation, ip):
        transformation.Transform(ip, self.point)
        p = self.point.GetDataArray()
        output.SetSize(self.dimension)
        values = output.GetDataArray()
        values[0] = self.scale * cos(kappa*p[0])*sin(kappa*p[1])
        values[1] = self.scale * cos(kappa*p[1])*sin(kappa*p[0])
        if self.dimension == 3:
            values[2] = 0.0


CALLBACK_COMPILE_SECONDS = perf_counter() - CALLBACK_COMPILE_STARTED


class PythonCoefficient(mfem.VectorPyCoefficient):
    def __init__(self, dimension, scale):
        super().__init__(dimension)
        self.scale = scale

    def EvalValue(self, p):
        values = [self.scale*cos(kappa*p[0])*sin(kappa*p[1]),
                  self.scale*cos(kappa*p[1])*sin(kappa*p[0])]
        if len(p) == 3:
            values.append(0.0)
        return np.array(values)


SETUP_SECONDS = perf_counter() - IMPORT_STARTED


def run(backend="bridge", ref_levels=None, order=1, save=False):
    coefficient = NumbaCoefficient if backend == "bridge" else PythonCoefficient
    set_bc, static_cond, hybridization = True, False, False
    meshfile = str(Path(__file__).resolve().parents[2] / "data" / "star.mesh")
    mesh = mfem.Mesh(meshfile, 1, 1)
    dim = mesh.Dimension()
    if ref_levels is None:
        ref_levels = int(np.floor(np.log(25000./mesh.GetNE())/np.log(2.)/dim))
    for x in range(ref_levels):
        mesh.UniformRefinement()

    fec = mfem.RT_FECollection(order-1, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec)

    print("Number of finite element unknows : " + str(fespace.GetTrueVSize()))

    ess_tdof_list = intArray()
    if mesh.bdr_attributes.Size():
        ess_bdr = intArray(mesh.bdr_attributes.Max())
        if set_bc:
            ess_bdr.Assign(1)
        else:
            ess_bdr.Assign(0)
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

    b = mfem.LinearForm(fespace)
    f = coefficient(dim, 1.0 + 2.0 * kappa * kappa)
    dd = mfem.VectorFEDomainLFIntegrator(f)
    b.AddDomainIntegrator(dd)
    start = perf_counter()
    b.Assemble()
    assembly_seconds = perf_counter() - start
    rhs = b.GetDataArray().copy()

    x = mfem.GridFunction(fespace)
    F = coefficient(dim, 1.0)
    start = perf_counter()
    x.ProjectCoefficient(F)
    projection_seconds = perf_counter() - start
    projected = x.GetDataArray().copy()

    alpha = mfem.ConstantCoefficient(1.0)
    beta = mfem.ConstantCoefficient(1.0)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DivDivIntegrator(alpha))
    a.AddDomainIntegrator(mfem.VectorFEMassIntegrator(beta))


    if (static_cond):
        a.EnableStaticCondensation()
    elif (hybridization):
        hfec = mfem.DG_Interface_FECollection(order-1, dim)
        hfes = mfem.FiniteElementSpace(mesh, hfec)
        a.EnableHybridization(hfes, mfem.NormalTraceJumpIntegrator(),
                              ess_tdof_list)
    a.Assemble()


    A = mfem.OperatorPtr()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
    # Here, original version calls hegith, which is not
    # defined in the header...!?
    print("Size of linear system: " + str(A.Height()))

    # 10. Solve
    AA = mfem.OperatorHandle2SparseMatrix(A)
    M = mfem.GSSmoother(AA)
    mfem.PCG(AA, M, B, X, 1, 10000, 1e-20, 0.0)

    # 11. Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, x)

    error = x.ComputeL2Error(F)
    print("|| F_h - F ||_{L^2} = " + str(error))

    if save:
        mesh.Print('refined.mesh', 16)
        x.Save('sol.gf', 16)
    result = dict(error=error, rhs=rhs, projected=projected,
                  solution=x.GetDataArray().copy(), assembly=assembly_seconds,
                  projection=projection_seconds)
    if backend == "bridge":
        f.close()
        F.close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("bridge", "python"), default="bridge")
    parser.add_argument("--refine", type=int, default=None)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    start = perf_counter()
    fresh = run(args.backend, args.refine, args.order, args.save)
    print(json.dumps(dict(phase="fresh", backend=args.backend,
                          elapsed=perf_counter()-start, setup_including_compilation=SETUP_SECONDS,
                          callback_compilation=CALLBACK_COMPILE_SECONDS,
                          assembly=fresh["assembly"], projection=fresh["projection"])))
    if args.compare:
        start = perf_counter()
        warmed = run(args.backend, args.refine, args.order)
        print(json.dumps(dict(phase="warmed", backend=args.backend,
                              elapsed=perf_counter()-start,
                              assembly=warmed["assembly"], projection=warmed["projection"])))
        start = perf_counter()
        reference = run("python", args.refine, args.order)
        print(json.dumps(dict(phase="reference", backend="python",
                              elapsed=perf_counter()-start,
                              assembly=reference["assembly"], projection=reference["projection"])))
        for field in ("rhs", "projected", "solution", "error"):
            np.testing.assert_allclose(fresh[field], reference[field], rtol=1e-10, atol=1e-11)
            np.testing.assert_allclose(warmed[field], reference[field], rtol=1e-10, atol=1e-11)
        print("Python and NSB fields, RHS, solution, and L2 error agree.")


if __name__ == "__main__":
    main()
