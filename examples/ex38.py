'''
   PyMFEM example 38

   See c++ version in the MFEM library for more detail

   Sample runs:
       python ex38.py
       python ex38.py -i volumetric1d
       python ex38.py -i surface2d
       python ex38.py -i surface2d -o 4 -r 5
       python ex38.py -i volumetric2d
       python ex38.py -i volumetric2d -o 4 -r 5
       python ex38.py -i surface3d
       python ex38.py -i surface3d -o 4 -r 5
       python ex38.py -i volumetric3d
       python ex38.py -i volumetric3d -o 4 -r 5

'''
import mfem.ser as mfem
from mfem.ser import intArray, doubleArray

import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, cos, array, pi, sqrt, floor, nan

import numpy as np
from numba import njit
from numba.types import float64


class IntegrationType:
    Volumetric1D = 0
    Surface2D = 1
    Volumetric2D = 2
    Surface3D = 3
    Volumetric3D = 4


visualization = True
itype = None


def make_lvlset_coeff():
    if itype == IntegrationType.Volumetric1D:
        @mfem.jit.scalar
        def func(X):
            return .55 - X[0]
    elif itype == IntegrationType.Surface2D:
        @mfem.jit.scalar
        def func(X):
            return 1. - (X[0]**2 + X[1]**2)
    elif itype == IntegrationType.Volumetric2D:
        @mfem.jit.scalar
        def func(X):
            return 1. - ((X[0]/1.5)**2 + (X[1]/.75)**2)
    elif itype == IntegrationType.Surface3D:
        @mfem.jit.scalar
        def func(X):
            return 1. - (X[0]**2 + X[1]**2 + X[2]**2)
    elif itype == IntegrationType.Volumetric3D:
        @mfem.jit.scalar
        def func(X):
            return 1. - ((X[0]/1.5)**2 + (X[1]/0.75)**2 + (X[2]/.5)**2)
    else:
        @njit(float64(float64[:]))
        def func(X):
            return 1.
    return func


def make_integrand_coeff():
    if itype == IntegrationType.Volumetric1D:
        @mfem.jit.scalar
        def func(X):
            return 1.
    elif itype == IntegrationType.Surface2D:
        @mfem.jit.scalar
        def func(X):
            return 3. * X[0]**2 - X[1]**2
    elif itype == IntegrationType.Volumetric2D:
        @mfem.jit.scalar
        def func(X):
            return 1.
    elif itype == IntegrationType.Surface3D:
        @mfem.jit.scalar
        def func(X):
            return 4. - 3. * X[0]**2 + 2. * X[1]**2 - X[2]**2
    elif itype == IntegrationType.Volumetric3D:
        @mfem.jit.scalar
        def func(X):
            return 1.
    else:
        @mfem.jit.scalar
        def func(X):
            return 0.
    return func


def Surface():
    if itype == IntegrationType.Volumetric1D:
        return 1.
    elif itype == IntegrationType.Surface2D:
        return 2. * pi
    elif itype == IntegrationType.Volumetric2D:
        return 7.26633616541076
    elif itype == IntegrationType.Surface3D:
        return 40. / 3. * pi
    elif itype == IntegrationType.Volumetric3D:
        return 9.90182151329315
    else:
        return 0.


def Volume():
    if itype == IntegrationType.Volumetric1D:
        return .55
    elif itype == IntegrationType.Surface2D:
        return nan
    elif itype == IntegrationType.Volumetric2D:
        return 9/8*pi
    elif itype == IntegrationType.Surface3D:
        return nan
    elif itype == IntegrationType.Volumetric3D:
        return 3/4*pi
    else:
        return 0.


class SIntegrationRule(mfem.PyIntegrationRule):
    def __init__(self, Order, LvlSet, lsOrder, mesh):
        super(SIntegrationRule, self).__init__()

        self.Weights = mfem.DenseMatrix()
        self.SurfaceWeights = mfem.DenseMatrix()
        self.dim = mesh.Dimension()

        MFIRs = mfem.MomentFittingIntRules(Order, LvlSet, lsOrder)
        Tr = mesh.GetElementTransformation(0)

        ir = mfem.IntegrationRule()
        MFIRs.GetSurfaceIntegrationRule(Tr, ir)
        if self.dim > 1:
            self.Weights.SetSize(ir.GetNPoints(), mesh.GetNE())
        else:
            self.Weights.SetSize(2, mesh.GetNE())

        self.SurfaceWeights.SetSize(ir.GetNPoints(), mesh.GetNE())

        w = mfem.Vector()
        MFIRs.GetSurfaceWeights(Tr, ir, w)
        self.SurfaceWeights.SetCol(0, w)

        self.SetSize(ir.GetNPoints())

        for ip in range(self.GetNPoints()):
            self.IntPoint(ip).index = ip
            intp = self.IntPoint(ip)
            intp.x = ir.IntPoint(ip).x
            intp.y = ir.IntPoint(ip).y
            intp.z = ir.IntPoint(ip).z
            if self.dim > 1:
                self.Weights[ip, 0] = ir.IntPoint(ip).weight
            else:
                self.Weights[0, 0] = ir.IntPoint(ip).x
                self.Weights[1, 0] = ir.IntPoint(ip).weight

        for elem in range(mesh.GetNE()):
            Tr = mesh.GetElementTransformation(elem)
            MFIRs.GetSurfaceIntegrationRule(Tr, ir)

            w = mfem.Vector()
            MFIRs.GetSurfaceWeights(Tr, ir, w)
            self.SurfaceWeights.SetCol(elem, w)

            for ip in range(self.GetNPoints()):
                if self.dim > 1:
                    self.Weights[ip, elem] = ir.IntPoint(ip).weight
                else:
                    self.Weights[0, elem] = ir.IntPoint(ip).x
                    self.Weights[1, elem] = ir.IntPoint(ip).weight

    def SetElementinclSurfaceWeight(self, Element):
        if self.dim == 1:
            intp = self.IntPoint(0)
            intp.x = self.Weights[0, Element]
            intp.weight = self.Weights[1, Element]
            print(str(intp.x) + " " + str(Element))
        else:
            for ip in range(self.GetNPoints()):
                intp = self.IntPoint(ip)
                intp.weight = self.Weights[ip, Element] * \
                    self.SurfaceWeights[ip, Element]

    def SetElement(self, Element):
        if self.dim == 1:
            intp = IntPoint(0)
            intp.x = self.Weights[0, Element]
            intp.weight = Weights[1, Element]
        else:
            for ip in range(self.GetNPoints()):
                intp = self.IntPoint(ip)
                intp.weight = self.Weights[ip, Element]


class CIntegrationRule(mfem.PyIntegrationRule):
    def __init__(self, Order, LvlSet, lsOrder, mesh):
        super(CIntegrationRule, self).__init__()
        self.dim = mesh.Dimension()
        self.Weights = mfem.DenseMatrix()

        MFIRs = mfem.MomentFittingIntRules(Order, LvlSet, lsOrder)

        Tr = mesh.GetElementTransformation(0)

        ir = mfem.IntegrationRule()
        MFIRs.GetVolumeIntegrationRule(Tr, ir)
        if self.dim > 1:
            self.Weights.SetSize(ir.GetNPoints(), mesh.GetNE())
        else:
            self.Weights.SetSize(2 * ir.GetNPoints(), mesh.GetNE())

        self.SetSize(ir.GetNPoints())
        for ip in range(self.GetNPoints()):
            self.IntPoint(ip).index = ip
            intp = self.IntPoint(ip)
            intp.x = ir.IntPoint(ip).x
            intp.y = ir.IntPoint(ip).y
            intp.z = ir.IntPoint(ip).z
            if self.dim > 1:
                self.Weights[ip, 0] = ir.IntPoint(ip).weight
            else:
                self.Weights[2 * ip, 0] = ir.IntPoint(ip).x
                self.Weights[2 * ip + 1, 0] = ir.IntPoint(ip).weight

        for elem in range(mesh.GetNE()):
            Tr = mesh.GetElementTransformation(elem)
            MFIRs.GetVolumeIntegrationRule(Tr, ir)
            for ip in range(self.GetNPoints()):
                if self.dim > 1:
                    self.Weights[ip, elem] = ir.IntPoint(ip).weight
                else:
                    self.Weights[2 * ip, elem] = ir.IntPoint(ip).x
                    self.Weights[2 * ip + 1, elem] = ir.IntPoint(ip).weight

    def SetElement(self, Element):
        if self.dim == 1:
            for ip in range(self.GetNPoints()):
                intp = self.IntPoint(ip)
                intp.x = self.Weights[2 * ip, Element]
                intp.weight = self.Weights[2 * ip + 1, Element]
        else:
            for ip in range(self.GetNPoints()):
                intp = self.IntPoint(ip)
                intp.weight = self.Weights[ip, Element]


class SurfaceLFIntegrator(mfem.PyLinearFormIntegrator):
    def __init__(self, q, levelset, ir):
        self.Q = q
        self.LevelSet = levelset
        self.SIntRule = ir
        self.shape = mfem.Vector()
        super(SurfaceLFIntegrator, self).__init__(ir)

    def AssembleRHSElementVect(self, el, Tr, elvect):
        dof = el.GetDof()
        self.shape.SetSize(dof)
        elvect.SetSize(dof)
        elvect.Assign(0.)

        # Update the surface integration rule for the current element
        self.SIntRule.SetElementinclSurfaceWeight(Tr.ElementNo)

        for ip in range(self.SIntRule.GetNPoints()):
            Tr.SetIntPoint(self.SIntRule.IntPoint(ip))
            val = Tr.Weight() * self.Q.Eval(Tr, self.SIntRule.IntPoint(ip))
            el.CalcShape(self.SIntRule.IntPoint(ip), self.shape)
            mfem.add_vector(elvect, self.SIntRule.IntPoint(
                ip).weight * val, self.shape, elvect)


class SubdomainLFIntegrator(mfem.PyLinearFormIntegrator):
    def __init__(self, q, levelset, ir):
        self.shape = mfem.Vector()
        self.CIntRule = ir
        self.LevelSet = levelset
        self.Q = q
        super(SubdomainLFIntegrator, self).__init__(ir)

    def AssembleRHSElementVect(self, el, Tr, elvect):
        dof = el.GetDof()
        self.shape.SetSize(dof)
        elvect.SetSize(dof)
        elvect.Assign(0.)

        # Update the subdomain integration rule
        self.CIntRule.SetElement(Tr.ElementNo)
        for ip in range(self.CIntRule.GetNPoints()):
            Tr.SetIntPoint(self.CIntRule.IntPoint(ip))
            val = Tr.Weight() * self.Q.Eval(Tr, self.CIntRule.IntPoint(ip))
            el.CalcPhysShape(Tr, self.shape)
            mfem.add_vector(elvect, self.CIntRule.IntPoint(
                ip).weight * val, self.shape, elvect)


def run(ref_levels=3, order=2):
    if itype == IntegrationType.Volumetric1D:
        mesh = mfem.Mesh("../data/inline-segment.mesh")
    elif (itype == IntegrationType.Surface2D or
          itype == IntegrationType.Volumetric2D):
        mesh = mfem.Mesh(2, 4, 1, 0, 2)
        mesh.AddVertex(-1.6, -1.6)
        mesh.AddVertex(1.6, -1.6)
        mesh.AddVertex(1.6, 1.6)
        mesh.AddVertex(-1.6, 1.6)
        mesh.AddQuad(0, 1, 2, 3)
        mesh.FinalizeQuadMesh(1, 0, True)
    elif (itype == IntegrationType.Surface3D or
          itype == IntegrationType.Volumetric3D):
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
    else:
        assert False, "Unknown integration type"

    for i in range(ref_levels):
        mesh.UniformRefinement()

    # 3. Define the necessary finite element space on the mesh.
    fe_col = mfem.H1_FECollection(1, mesh.Dimension())
    fespace = mfem.FiniteElementSpace(mesh, fe_col)

    # 4. Construction Coefficients for the level set and the integrand.
    levelset = make_lvlset_coeff()
    u = make_integrand_coeff()

    # 5. Define the necessary Integration rules on element 0.
    Tr = mesh.GetElementTransformation(0)

    sir = SIntegrationRule(order, levelset, 2, mesh)
    if (itype == IntegrationType.Volumetric1D or
        itype == IntegrationType.Volumetric2D or
            itype == IntegrationType.Volumetric3D):
        cir = CIntegrationRule(order, levelset, 2, mesh)

    # 6. Define and assemble the linear forms on the finite element space.
    surface = mfem.LinearForm(fespace)
    volume = mfem.LinearForm(fespace)

    surface.AddDomainIntegrator(SurfaceLFIntegrator(u, levelset, sir))
    surface.Assemble()

    if (itype == IntegrationType.Volumetric1D or
        itype == IntegrationType.Volumetric2D or
            itype == IntegrationType.Volumetric3D):
        volume.AddDomainIntegrator(SubdomainLFIntegrator(u, levelset, cir))
        volume.Assemble()

    # 7. Print information, computed values and errors to the console.
    qorder = 0
    nbasis = 2 * (order + 1) + (order * (order + 1) // 2)

    irs = mfem.IntegrationRules(0, mfem.Quadrature1D.GaussLegendre)
    ir = irs.Get(mfem.Geometry.SQUARE, qorder)

    while (ir.GetNPoints() <= nbasis):
        ir = irs.Get(mfem.Geometry.SQUARE, qorder)
        qorder += 1

    print("============================================")
    if itype != IntegrationType.Volumetric1D:
        print("Mesh size dx:                       " +
              "{:g}".format(3.2 / 2**ref_levels))
    else:
        print("Mesh size dx:                       " +
              "{:g}".format(25 / 2**ref_levels))

    if (itype == IntegrationType.Surface2D or
            itype == IntegrationType.Volumetric2D):
        print("Number of div free basis functions: " + str(nbasis))
        print("Number of quadrature points:        " + str(ir.GetNPoints()))

    print("============================================")
    print("Computed value of surface integral: " +
          "{:.2e}".format(surface.Sum()))
    print("True value of surface integral:     " + "{:.2e}".format(Surface()))

    print("Absolut Error (Surface):            " +
          "{:.2e}".format(abs(surface.Sum() - Surface())))
    print("Relative Error (Surface):           " +
          "{:.2e}".format(abs(surface.Sum() - Surface()) / Surface()))

    if (itype == IntegrationType.Volumetric1D or
        itype == IntegrationType.Volumetric2D or
            itype == IntegrationType.Volumetric3D):

        print("--------------------------------------------")
        print("Computed value of volume integral:  " +
              "{:.2e}".format(volume.Sum()))
        print("True value of volume integral:      " +
              "{:.2e}".format(Volume()))
        print("Absolut Error (Volume):             " +
              "{:.2e}".format(abs(volume.Sum() - Volume())))
        print("Relative Error (Volume):            " +
              "{:.2e}".format(abs(volume.Sum() - Volume()) / Volume()))

    print("============================================")

    # 8. Plot the level-set function on a high order finite element space.
    if visualization:
        fe_coll2 = mfem.H1_FECollection(5, mesh.Dimension())
        fespace2 = mfem.FiniteElementSpace(mesh, fe_coll2)
        levelset_coeff = make_lvlset_coeff()

        lgf = mfem.GridFunction(fespace2)
        lgf.ProjectCoefficient(levelset_coeff)

        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock << "solution\n" << mesh << lgf
        sol_sock.flush()
        sol_sock << "keys pppppppppppppppppppppppppppcmmlRj\n"
        sol_sock << "levellines " << "0." << " " << "0." << " " << 1 << "\n"
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(
        description='Ex38 (Cut-Volume and Cut-Surface Integration)')

    parser.add_argument("-o", "--order",
                        action='store', default=2, type=int,
                        help="Order (degree) of the finite elements.")
    parser.add_argument("-r", "--refine",
                        action='store', default=3, type=int,
                        help="Number of meh refinements")
    parser.add_argument("-i", "--integrationtype",
                        action='store', default="surface2d", type=str,
                        help="IntegrationType to demonstrate")
    parser.add_argument("-no-vis", "--no-visualization",
                        action='store_true', default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    globals()["visualization"] = not args.no_visualization

    if args.integrationtype.capitalize().startswith("Volumetric1"):
        globals()["itype"] = IntegrationType.Volumetric1D
    elif args.integrationtype.capitalize().startswith("Volumetric2"):
        globals()["itype"] = IntegrationType.Volumetric2D
    elif args.integrationtype.capitalize().startswith("Volumetric3"):
        globals()["itype"] = IntegrationType.Volumetric3D
    elif args.integrationtype.capitalize().startswith("Surface2"):
        globals()["itype"] = IntegrationType.Surface2D
    elif args.integrationtype.capitalize().startswith("Surface3"):
        globals()["itype"] = IntegrationType.Surface3D
    else:
        assert False, "Unknown integration type"

    if not hasattr(mfem, "MomentFittingIntRules"):
        assert False, "MFEM is not built with Lapack"

    run(ref_levels=args.refine,
        order=args.order)
