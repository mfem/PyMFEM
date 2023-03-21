'''
   -------------------------------------------------------------
   Field Interp Miniapp: Transfer a grid function between meshes
   -------------------------------------------------------------

   See C++ example for documentation

   Sample runs:
      python field_interp.py -v
      python field_interp.py -fts 3 -ft 0 -v
      python field_interp.py -m1 triple-pt-1.mesh -s1 triple-pt-1.gf -m2 triple-pt-2.mesh -ft 1 -v
      python field_interp.py -m2 ../../data/amr-quad-q2.mesh -ft 0 -r 1 -v
'''

import os
from os.path import expanduser, join
import numpy as np

try:
    import mpi4py
    import mfem.par as mfem
except:
    import mfem.ser as mfem

def run(order=3,
        meshfile1='',
        meshfile2='',
        solution1='',
        visualization=True,
        src_fieldtype=0,
        src_ncomp=1,
        ref_levels=0,
        fieldtype=-1):

    # Input meshes.
    mesh_1 = mfem.Mesh(meshfile1, 1, 1, False)
    mesh_2 = mfem.Mesh(meshfile2, 1, 1, False)

    dim = mesh_1.Dimension()

    errtxt = "Source and target meshes must be in the same dimension."
    assert dim == mesh_2.Dimension(), errtxt
    assert dim <= 2, "GSLIB requires a 2D or a 3D mesh"

    for i in range(ref_levels):
        mesh_2.UniformRefinement()

    if mesh_1.GetNodes() is None:
        mesh_1.SetCurvature(1)
    if mesh_2.GetNodes() is None:
        mesh_2.SetCurvature(1)

    mesh_poly_deg = mesh_2.GetNodes().FESpace().GetElementOrder(0)

    print("Source mesh curvature: " + mesh_1.GetNodes().OwnFEC().Name())
    print("Target mesh curvature: " + mesh_2.GetNodes().OwnFEC().Name())

    if src_fieldtype < 0:
        func_source = mfem.GridFunction(mesh_1, solution1)
        src_vdim = func_source.FESpace().GetVDim()
    elif src_fieldtype == 0:
        src_fec = mfem.H1_FECollection(order, dim)
        src_vdim = src_ncomp
    elif src_fieldtype == 1:
        src_fec = mfem.L2_FECollection(order, dim)
        src_vdim = src_ncomp
    elif src_fieldtype == 2:
        src_fec = mfem.RT_FECollection(order, dim)
        src_ncomp = 1
        src_vdim = dim
    elif src_fieldtype == 3:
        src_fec = mfem.ND_FECollection(order, dim)
        src_ncomp = 1
        src_vdim = dim
    else:
        assert False, "Invalid FECollection type."

    sdim = mesh_1.SpaceDimension()

    @mfem.jit.scalar()
    def scalar_func(x):
        res = 0.0
        for d in range(len(x)):
            res += x[d] * x[d]
        return res

    @mfem.jit.vector(vdim=src_vdim)
    def vector_func(x):
        F = np.zeros(vdim, dtype=float)
        for d in range(len(x)):
            F[0] += x[d] * x[d]
        for i in range(1, src_vdim):
            F[i] = (i + 1) * (-1)**i * F[0]
        return F

    if src_fieldtype > -1:
        src_fes = mfem.FiniteElementSpace(mesh_1, src_fec, src_ncomp)
        func_source = mfem.GridFunction(src_fes)

        # Project the grid function using numba compiled
        # VectorFunctionCoefficient.
        func_source.ProjectCoefficient(vector_func)

    if visualization:
        sout1 = mfem.socketstream("localhost", 19916)
        if not sout1.good():
            print("Unable to connect to GLVis server at localhost 19916")
        else:
            sout1 << "solution\n" << mesh_1 << func_source
            sout1 << "window_title 'Source mesh and solution'"
            sout1 << "window_geometry 0 0 600 600"
            if dim == 2:
                sout1 << "keys RmjAc"
            if dim == 3:
                sout1 << "keys mA\n"
            sout1.flush()

    gt = mesh_2.GetNodalFESpace().GetFE(0).GetGeomType()
    assert gt != mfem.Geometry.PRISM, "Wedge elements are not currently supported."
    assert mesh_2.GetNumGeometries(
        mesh_2.Dimension()) == 1, "Mixed meshes are not currently supported."

    # Ensure the source grid function can be transferred using GSLIB-FindPoints.
    fec_in = func_source.FESpace().FEColl()
    print("Source FE collection: " + fec_in.Name())

    if src_fieldtype < 0:
        if fec_in.Name().startswith('H1'):
            src_fieldtype = 0
        elif fec_in.Name().startswith('L2'):
            src_fieldtype = 1
        elif fec_in.Name().startswith('RT'):
            src_fieldtype = 2
        elif fec_in.Name().startswith('ND'):
            src_fieldtype = 3
        else:
            assert False, "GridFunction type not supported yet."

    if fieldtype < 0:
        fieldtype = src_fieldtype

    tar_vdim = src_vdim
    if fieldtype == 0:
        tar_fec = mfem.H1_FECollection(order, dim)
        tar_vdim = dim if src_fieldtype > 1 else src_vdim
    elif fieldtype == 1:
        tar_fec = mfem.L2_FECollection(order, dim)
        tar_vdim = dim if src_fieldtype > 1 else src_vdim
    elif fieldtype == 2:
        tar_fec = mfem.RT_FECollection(order, dim)
        tar_vdim = 1
        assert src_fieldtype != 1, "Cannot interpolate a scalar grid function to a vector"
    elif fieldtype == 3:
        tar_fec = mfem.ND_FECollection(order, dim)
        tar_vdim = 1
        assert src_fieldtype != 1, "Cannot interpolate a scalar grid function to a vector"
    else:
        assert False, "GridFunction type not supported."

    print("Target FE collection: " + tar_fec.Name())
    tar_fes = mfem.FiniteElementSpace(mesh_2, tar_fec, tar_vdim)
    func_target = mfem.GridFunction(tar_fes)

    NE = mesh_2.GetNE()
    nsp = tar_fes.GetFE(0).GetNodes().GetNPoints()
    tar_ncomp = func_target.VectorDim()

    # Generate list of points where the grid function will be evaluated.

    if fieldtype == 0 and order == mesh_poly_deg:
        vxyz = mesh_2.GetNodes()
    else:
        vxyz0 = np.empty(nsp * NE * dim)

        for i in range(NE):
            fe = tar_fes.GetFE(i)
            ir = fe.GetNodes()
            et = tar_fes.GetElementTransformation(i)

            pos = mfem.DenseMatrix()
            et.Transform(ir, pos)
            pos0 = pos.GetDataArray()

            vxyz0[i * nsp:i * nsp + nsp] = pos0[0, :]
            vxyz0[i * nsp + NE * nsp:i * nsp + NE * nsp + nsp] = pos0[1, :]
            if dim == 3:
                vxyz0[i * nsp + 2 * NE * nsp: i * nsp +
                      2 * NE * nsp + nsp] = pos0[2, :]

        vxyz = mfem.Vector(vxyz0)

    nodes_cnt = vxyz.Size() // dim

    # Evaluate source grid function
    interp_vals = mfem.Vector(nodes_cnt * tar_ncomp)
    finder = mfem.FindPointsGSLIB()
    finder.Setup(mesh_1)
    finder.Interpolate(vxyz, func_source, interp_vals)

    interp_vals0 = interp_vals.GetDataArray()

    # Project the interpolated values to the target FiniteElementSpace.
    if fieldtype <= 1:
        if (fieldtype == 0 and order == mesh_poly_deg) or fieldtype == 1:
            func_target.Assign(interp_vals)
        else:  # H1 - but mesh order != GridFunction order
            elem_dof_vals = np.empty(nsp * tar_ncomp, dtype=float)
            func_target0 = func_target.GetDataArray()

            for i in range(mesh_2.GetNE()):
                vdofs = tar_fes.GetElementVDofs(i)
                for j in range(nsp):
                    for d in range(tar_ncomp):
                        elem_dof_vals[j + d * nsp] = interp_vals0[d *
                                                                  nsp * NE + i * nsp + j]

                func_target0[vdofs] = elem_dof_vals
    else:  # H(div) or H(curl)
        vals = mfem.Vector()
        elem_dof_vals = np.zeros(nsp * tar_ncomp, dtype=float)
        for i in range(mesh_2.GetNE()):
            vdofs = tar_fes.GetElementVDofs(i)
            vals.SetSize(len(vdofs))
            for j in range(nsp):
                for d in range(tar_ncomp):
                    # Arrange values byVDim
                    elem_dof_vals[j * tar_ncomp +
                                  d] = interp_vals0[d * nsp * NE + i * nsp + j]
            tar_fes.GetFE(i).ProjectFromNodes(elem_dof_vals,
                                              tar_fes.GetElementTransformation(
                                                  i),
                                              vals)
            func_target0[vdofs] = vals.GetDataArray()

    if visualization:
        sout1 = mfem.socketstream("localhost", 19916)
        if not sout1.good():
            print("Unable to connect to GLVis server at localhost 19916")
        else:
            sout1 << "solution\n" << mesh_2 << func_target
            sout1 << "window_title 'Target mesh and solution'"
            sout1 << "window_geometry 0 0 600 600"
            if dim == 2:
                sout1 << "keys RmjAc"
            if dim == 3:
                sout1 << "keys mA\n"
            sout1.flush()


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Field Interp')
    parser.add_argument('-m1', '--mesh1',
                        default="../../data/square01.mesh",
                        action='store', type=str,
                        help="Mesh file for the starting solution.")
    parser.add_argument('-m2', '--mesh2',
                        default="../../data/inline-tri.mesh",
                        action='store', type=str,
                        help="Mesh file for interpolation.")

    txt = ["(optional) GridFunction file compatible with src_mesh_file.",
           "Set src_fieldtype to -1 if this option is used."]
    parser.add_argument("-s1", "--solution1",
                        type=str, action='store',
                        default="must_be_provided_by_the_user.gf",
                        help="\n".join(txt),)

    txt = ["Source GridFunction type:",
           "0 - H1 (default), 1 - L2, 2 - H(div), 3 - H(curl)."]
    parser.add_argument("-fts", "--field-type-src",
                        default=0, type=int,
                        help="\n".join(txt))

    parser.add_argument("-nc", "--ncomp",
                        default=1, type=int,
                        help="Number of components for H1 or L2 GridFunctions.")
    parser.add_argument("-r", "--refine",
                        default=0, type=int,
                        help="Number of refinements of the interpolation mesh.")

    txt = ["Target GridFunction type: -1 - source GridFunction type (default),",
           "0 - H1, 1 - L2, 2 - H(div), 3 - H(curl)."]
    parser.add_argument("-ft", "--field-type",
                        default=-1, type=int,
                        help="\n".join(txt))

    parser.add_argument('-o', '--order',
                        action='store', default=3, type=int,
                        help="Order of the interpolated solution.")
    parser.add_argument('-vis', '--visualization',
                        action='store_true',
                        default=False,
                        help='Enable GLVis visualization')

    args = parser.parse_args()
    parser.print_options(args)

    order = args.order

    meshfile1 = expanduser(
        join(os.path.dirname(__file__), args.mesh1))
    meshfile2 = expanduser(
        join(os.path.dirname(__file__), args.mesh2))

    run(order=order,
        meshfile1=meshfile1,
        meshfile2=meshfile2,
        solution1=args.solution1,
        visualization=args.visualization,
        src_fieldtype=args.field_type_src,
        src_ncomp=args.ncomp,
        ref_levels=args.refine,
        fieldtype=args.field_type)
