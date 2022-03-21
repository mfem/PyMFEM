'''
   MFEM example 30

   See c++ version in the MFEM library for more detail 
'''
import mfem.ser as mfem
from mfem.ser import intArray
import os
from os.path import expanduser, join
import numpy as np
from numpy import sin, array

try:
    from numba import jit
except:
    assert False, "This example requires numba"


@mfem.jit.scalar()
def affine_coeff(p):
    if p[0] < 0.0:
        return 1 + p[0] + p[1]
    else:
        return 1


@mfem.jit.scalar()
def jump_coeff(p):
    norm = np.sqrt(p[0]**2 + p[1]**2)
    if norm > 0.4 and norm < 0.6:
        return 1.0
    else:
        return 5.0


@mfem.jit.scalar()
def singular_coeff(p):
    x = p[0]
    y = p[1]
    alpha = 1000.0
    xc = 0.75
    yc = 0.5
    r0 = 0.7
    r = np.sqrt((x - xc)**2.0 + (y - yc)**2)
    num = - (alpha - alpha**3 * (r**2 - r0**2))
    denom = (r * (alpha**2 * r0**2 + alpha**2 * r**2 -
                 2 * alpha**2 * r0 * r + 1.0))**2
    denom = max([denom, 1e-8])
    return num / denom


def run(nc_limit=1,
        order=1,
        meshfile='',
        visualization=True,
        max_elems=100000,
        enriched_order=5,
        osc_threshold=1e-3):

    mesh = mfem.Mesh(meshfile, 1, 1)

    if mesh.NURBSext:
        for i in range(2):
            mesh.UniformRefinement()
        mesh.SetCurvature(2)

    coeffrefiner = mfem.CoefficientRefiner(affine_coeff, order)

    if (visualization):
        sol_sock = mfem.socketstream("localhost", 19916)

    irs = [mfem.IntRules.Get(i, 2*order + enriched_order)
           for i in range(mfem.Geometry.NumGeom)]

    coeffrefiner.SetIntRule(irs)
    coeffrefiner.SetMaxElements(max_elems)
    coeffrefiner.SetThreshold(osc_threshold)
    coeffrefiner.SetNCLimit(nc_limit)
    coeffrefiner.PrintWarnings()

    coeffrefiner.PreprocessMesh(mesh)

    print()
    print("Function 0 (affine)")
    print("Number of Elements " + str(mesh.GetNE()))
    print("Osc error " + "{:g}".format(coeffrefiner.GetOsc()))

    # 8. Preprocess mesh to control osc (jump function).
    coeffrefiner.ResetCoefficient(jump_coeff)
    coeffrefiner.PreprocessMesh(mesh)

    print("")
    print("Function 1 (discontinuous)")
    print("Number of Elements " + str(mesh.GetNE()))
    print("Osc error " + "{:g}".format(coeffrefiner.GetOsc()))

    # 9. Preprocess mesh to control osc (singular function).
    coeffrefiner.ResetCoefficient(singular_coeff)
    coeffrefiner.PreprocessMesh(mesh)

    print("")
    print("Function 2 (singular)")
    print("Number of Elements " + str(mesh.GetNE()))
    print("Osc error " + "{:g}".format(coeffrefiner.GetOsc()))

    if visualization:
        sol_sock.precision(8)
        sol_sock << "mesh\n" << mesh
        sol_sock.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex30 (Coefficient Refiner)')
    parser.add_argument('-m', '--mesh',
                        default="star.mesh",
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument("-l", "--nc_limit",
                        action='store',
                        type=int,
                        default=1,
                        help="Maximum level of hanging nodes.")
    parser.add_argument('-me', '--max-elems',
                        action='store',
                        type=int, default=100000,
                        help="Stop after reaching this many elements",)
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        help='Enable GLVis visualization')
    parser.add_argument("-e", "--error",
                        action='store',
                        default=1e-3,
                        type=float,
                        help="relative data oscillation threshold.")
    parser.add_argument("-eo", "--enriched-order",
                        type=int,
                        default=5,
                        help="Enriched quadrature order.")

    args = parser.parse_args()
    parser.print_options(args)

    meshfile = expanduser(
        join(os.path.dirname(__file__), '..', 'data', args.mesh))
    print(meshfile)

    run(nc_limit=args.nc_limit,
        order=args.order,
        meshfile=meshfile,
        visualization=args.visualization,
        max_elems=args.max_elems,
        enriched_order=args.enriched_order,
        osc_threshold=args.error)
