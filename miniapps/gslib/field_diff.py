'''
   -------------------------------------------------------------
   Field Interp Miniapp: Compare grid functions on different meshes
   -------------------------------------------------------------

   See C++ example for documentation: https://github.com/mfem/mfem/blob/master/miniapps/gslib/field-diff.cpp

   Sample runs:
      python field_diff.py -vis
      python field_diff.py -m1 triple-pt-1.mesh -s1 triple-pt-1.gf -m2 triple-pt-2.mesh -s2 triple-pt-2.gf -p 200
'''

import os
from os.path import expanduser, join

import mfem.ser as mfem
import numpy as np


def field_diff(meshfile1='triple-pt-1.mesh',
               meshfile2='triple-pt-2.mesh',
               solution1='triple-pt-1.gf',
               solution2='triple-pt-2.gf',
               visualization=True,
               pts_cnt_1D=100
               ):
    """
    Compare two grid functions on two different meshes.
    Args:
        meshfile1: The mesh file for the first grid function.
        meshfile2: The mesh file for the second grid function.
        solution1: The grid function for the first mesh.
        solution2: The grid function for the second mesh.
        visualization: Whether to visualize the meshes and the grid functions.
            Requires an active GLVis server on port 19916.
        pts_cnt_1D: Number of comparison points in one direction. Will use pts_cnt_1D^dim points in total.

    Returns:

    """
    # ll 69-88 of field-diff.cpp
    # Input meshes.
    mesh_1 = mfem.Mesh(meshfile1, 1, 1, False)
    mesh_2 = mfem.Mesh(meshfile2, 1, 1, False)
    dim = mesh_1.Dimension()
    assert dim == mesh_2.Dimension(), "Source and target meshes must be in the same dimension."
    assert dim > 1, "GSLIB requires a 2D or a 3D mesh"

    # The Nodes GridFunctions for each mesh are required.
    if mesh_1.GetNodes() is None:
        mesh_1.SetCurvature(1)
    if mesh_2.GetNodes() is None:
        mesh_2.SetCurvature(1)
    mesh_poly_deg = mesh_1.GetNodes().FESpace().GetOrder(0)
    print(f"Mesh curvature: {mesh_1.GetNodes().OwnFEC().Name()} of order {mesh_poly_deg}")

    #  The visualization of the difference assumes byNODES ordering.
    if mesh_1.GetNodes().FESpace().GetOrdering() == mfem.Ordering.byVDIM:
        mesh_1.SetCurvature(mesh_poly_deg, False, dim, mfem.Ordering.byNODES)
    if mesh_2.GetNodes().FESpace().GetOrdering() == mfem.Ordering.byVDIM:
        mesh_2.SetCurvature(mesh_poly_deg, False, dim, mfem.Ordering.byNODES)

    # ll 88-98 of field-diff.cpp
    # Mesh bounding box
    assert mesh_poly_deg > 0, "The order of the mesh must be positive."
    pos_min, pos_max = mesh_1.GetBoundingBox(mesh_poly_deg)
    print(f"Generating equidistant points for:")
    print(f"  x in [{pos_min[0]}, {pos_max[0]}]")
    print(f"  y in [{pos_min[1]}, {pos_max[1]}]")
    if dim == 3:
        print(f"  z in [{pos_min[2]}, {pos_max[2]}]")

    # ll 99-103 of field-diff.cpp
    # Initialize grid functions using the provided solution files
    func_1 = mfem.GridFunction(mesh_1, solution1)
    func_2 = mfem.GridFunction(mesh_2, solution2)

    # ll 104-135 of field-diff.cpp
    if visualization:
        _visualize_mesh(func_1, mesh_1, dim, title="Source mesh and solution")
        _visualize_mesh(func_2, mesh_2, dim, title="Target mesh and solution")

    # ll 136-164 of field-diff.cpp. Note that this could be simplified in Python using e.g., np.linspace and np.meshgrid
    # Generate equidistant points in physical coordinates over the whole mesh.
    # Note that some points might be outside, if the mesh is not a box.
    # Note also that all tasks search the same points (not mandatory).
    pts_cnt = pts_cnt_1D ** dim
    vxyz = np.zeros(pts_cnt * dim)
    element_type = mfem.L2_QuadrilateralElement if dim == 2 else mfem.L2_HexahedronElement
    el = element_type(pts_cnt_1D - 1, mfem.BasisType.ClosedUniform)
    ir = el.GetNodes()
    for i in range(ir.GetNPoints()):
        ip = ir.IntPoint(i)
        vxyz[i] = pos_min[0] + ip.x * (pos_max[0] - pos_min[0])
        vxyz[pts_cnt + i] = pos_min[1] + ip.y * (pos_max[1] - pos_min[1])
        if dim == 3:
            vxyz[2 * pts_cnt + i] = pos_min[2] + ip.z * (pos_max[2] - pos_min[2])

    # ll 165-173 of field-diff.cpp
    interpolation_points = mfem.Vector(vxyz)
    interp_vals_1 = _interpolate_grid_function(func_1, mesh_1, pts_cnt, interpolation_points)  # first solution
    interp_vals_2 = _interpolate_grid_function(func_2, mesh_2, pts_cnt, interpolation_points)  # second solution

    # ll 175-208 of field-diff.cpp
    # Compute differences between the two sets of values
    diffs = np.abs(interp_vals_1.GetDataArray() - interp_vals_2.GetDataArray())
    avg_diff = np.mean(diffs)
    max_diff = np.max(diffs)

    # Compute the average distance between the two meshes if they have the same number of nodes
    n1 = mesh_1.GetNodes()
    n2 = mesh_2.GetNodes()
    if n1.Size() == n2.Size():
        avg_dist = 0.0
        node_cnt = n1.Size() / dim
        nd1 = n1.GetData()
        nd2 = n2.GetData()
        for i in range(node_cnt):
            diff_i = 0.0
            for d in range(dim):
                j = i + d * node_cnt
                diff_i += (nd1[j] - nd2[j]) * (nd1[j] - nd2[j])
            avg_dist += np.sqrt(diff_i)
            # np.linalg.norm variant:
            # avg_dist += np.linalg.norm(nd1[i:i+dim] - nd2[i:i+dim])
        avg_dist /= node_cnt
    else:
        avg_dist = -1.0

    print(f"Avg position difference: {avg_dist}\n"
          f"Searched {pts_cnt} points.\n"
          f"Max diff: {max_diff}\n"
          f"Avg diff: {avg_diff}")

    # ll 210-223 of field-diff.cpp
    # Interpolate the first solution onto the second mesh and compute the difference
    m1_fes = mfem.FiniteElementSpace(mesh_1, mesh_1.GetNodes().FESpace().FEColl())
    diff = mfem.GridFunction(m1_fes)
    f1c = mfem.GridFunctionCoefficient(func_1)
    diff.ProjectDiscCoefficient(f1c, mfem.GridFunction.ARITHMETIC)
    vxyz = mesh_1.GetNodes()
    nodes_cnt = int(vxyz.Size() / dim)
    interp_vals_2 = _interpolate_grid_function(func=func_2, mesh=mesh_2, pts_cnt=nodes_cnt, interpolation_points=vxyz)
    for n in range(nodes_cnt):
        diff[n] = np.abs(diff[n] - interp_vals_2[n])

    # ll 225-238 of field-diff.cpp
    if visualization:
        _visualize_mesh(diff, mesh_1, dim, title="Difference",
                        geometry="1200 0 600 600", keys="RmjAcpppppppppppppppppppppp")

    # ll 240-250 of field-diff.cpp
    coeff1 = mfem.ConstantCoefficient(1.0)
    lf_integ = mfem.DomainLFIntegrator(coeff1)
    lf = mfem.LinearForm(func_1.FESpace())
    lf.AddDomainIntegrator(lf_integ)
    lf.Assemble()
    vol_diff = np.dot(diff.GetDataArray(), lf.GetDataArray())
    print("Vol diff:", vol_diff)


def _interpolate_grid_function(func, mesh, pts_cnt: int, interpolation_points: mfem.Vector):
    # Evaluate source grid function
    interp_vals = mfem.Vector(pts_cnt)
    finder = mfem.FindPointsGSLIB()
    finder.Setup(mesh)
    finder.Interpolate(interpolation_points, func, interp_vals)
    return interp_vals


def _visualize_mesh(func, mesh, dim, title: str = "Mesh and solution",
                    vishost: str = "localhost", visport: int = 19916,
                    geometry: str = "0 0 600 600", keys: str = "RmjAc"):
    sout = mfem.socketstream(vishost, visport)
    if not sout:
        print(f"Unable to connect to GLVis server at {vishost}:{visport}")
        return
    # plot the source mesh and solution
    sout.precision(8)
    sout << "solution\n" << mesh << func
    sout << f"window_title '{title}'"
    sout << f"window_geometry {geometry}"
    if dim == 2:
        sout << f"keys {keys}"
    elif dim == 3:
        sout << "keys mA\n"
    sout.flush()


def main():
    # lines 36-67 of field-diff.cpp
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Field Interp')
    parser.add_argument('-m1', '--mesh1',
                        default="triple-pt-1.mesh",
                        action='store', type=str,
                        help="Mesh file for solution 1.")
    parser.add_argument('-m2', '--mesh2',
                        default="triple-pt-2.mesh",
                        action='store', type=str,
                        help="Mesh file for solution 2.")
    parser.add_argument("-s1", "--solution1",
                        type=str, action='store',
                        default="triple-pt-1.gf",
                        help="Grid function for solution 2.", )
    parser.add_argument("-s2", "--solution2",
                        type=str, action='store',
                        default="triple-pt-2.gf",
                        help="Grid function for solution 2.", )
    parser.add_argument('-p', "--points1D",
                        default=100, type=int,
                        help="Number of comparison points in one direction")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    meshfile1 = expanduser(join(os.path.dirname(__file__), args.mesh1))
    meshfile2 = expanduser(join(os.path.dirname(__file__), args.mesh2))
    solution1 = expanduser(join(os.path.dirname(__file__), args.solution1))
    solution2 = expanduser(join(os.path.dirname(__file__), args.solution2))
    points_1d = args.points1D
    visualization = args.visualization
    field_diff(meshfile1=meshfile1,
               meshfile2=meshfile2,
               solution1=solution1,
               solution2=solution2,
               pts_cnt_1D=points_1d,
               visualization=visualization)


if __name__ == '__main__':
    main()
