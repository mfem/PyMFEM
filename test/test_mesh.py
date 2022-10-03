from __future__ import print_function
import os
import sys
import numpy as np


def run_test1(mfem):
    print("Test mesh module")
    Nvert = 6
    Nelem = 8
    Nbelem = 2

    mesh = mfem.Mesh(2, Nvert, Nelem, 2, 3)
    tri_v = [[1.,  0.,  0.], [0.,  1.,  0.], [-1.,  0.,  0.],
             [0., -1.,  0.], [0.,  0.,  1.], [0.,  0., -1.]]
    tri_e = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4],
             [1, 0, 5], [2, 1, 5], [3, 2, 5], [0, 3, 5]]
    tri_l = [[1, 4], [1, 2]]

    for j in range(Nvert):
        mesh.AddVertex(tri_v[j])
    for j in range(Nelem):
        mesh.AddTriangle(tri_e[j], 1)
    for j in range(Nbelem):
        mesh.AddBdrSegment(tri_l[j], 1)

    mesh.FinalizeTriMesh(1, 1, True)

    print(mesh.GetEdgeVertices(1))
    print(mesh.GetFaceElements(1))

    assert mesh.GetBdrElementAdjacentElement(1) == (
        1, 0), "GetBdrElementadjacentElement failed"

    i = 0
    Tr = mesh.GetElementTransformation(i)
    print(Tr.Weight())
    Tr = mesh.GetBdrElementTransformation(i)
    print(Tr.Weight())
    Tr = mesh.GetFaceTransformation(i)
    print(Tr.Weight())
    Tr = mesh.GetEdgeTransformation(i)
    print(Tr.Weight())


def run_test2(mfem):
    import tracemalloc
    import linecache
    mesh = mfem.Mesh(6, 6, 6, "TETRAHEDRON")

    r = mesh.CartesianPartitioning([2, 2, 2], return_list=True)
    print(r)

    r = mesh.CartesianPartitioning([2, 2, 2])
    print(r)
    mesh.PrintWithPartitioning(r, "partiotioned.mesh")
    '''
    tracemalloc.start()
    def display_top(snapshot, key_type='lineno', limit=3):
        snapshot = snapshot.filter_traces((
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        ))
        top_stats = snapshot.statistics(key_type)

        print("Top %s lines" % limit)
        for index, stat in enumerate(top_stats[:limit], 1):
            frame = stat.traceback[0]
            # replace "/path/to/module/file.py" with "module/file.py"
            filename = os.sep.join(frame.filename.split(os.sep)[-2:])
            print("#%s: %s:%s: %.1f KiB"
                  % (index, filename, frame.lineno, stat.size / 1024))
            line = linecache.getline(frame.filename, frame.lineno).strip()
            if line:
                print('    %s' % line)

        other = top_stats[limit:]
        if other:
            size = sum(stat.size for stat in other)
            print("%s other: %.1f KiB" % (len(other), size / 1024))
        total = sum(stat.size for stat in top_stats)
        print("Total allocated size: %.1f KiB" % (total / 1024))

    snapshot = tracemalloc.take_snapshot()
    display_top(snapshot)
    '''
    '''
    d = mfem.intArray([2, 2, 2])
    dd = d.GetData()

    r = mesh.CartesianPartitioning(dd)
    result = mfem.intArray()
    result.MakeRef(r, mesh.GetNE())
    result.MakeDataOwner()
    '''

    '''    

    mesh.Print("sample.mesh")
    r = mesh.CartesianPartitioning([2, 2, 2])
    print(r)
    '''


def run_test3(mfem):
    m1 = mfem.Mesh_MakeCartesian3D(3, 3, 3, mfem.Element.TETRAHEDRON)
    m2 = mfem.Mesh_MakeCartesian3D(3, 3, 3, mfem.Element.TETRAHEDRON)

    for i in range(m2.GetNV()):
        vv = m2.GetVertexArray(i)
        vv += 2.5
    print(np.vstack(m1.GetVertexArray()))
    print(np.vstack(m2.GetVertexArray()))

    print("test if it fails (1)")
    try:
        m3 = mfem.Mesh_MakeMerged("error")
    except:
        import traceback
        print("exception happens")
        traceback.print_exc()
    print("test if it fails (2)")
    try:
        m3 = mfem.Mesh_MakeMerged(('hoge', m2))
    except:
        import traceback
        print("exception happens")
        traceback.print_exc()

    print("test if this one is ok (3)")
    m3 = mfem.Mesh_MakeMerged((m1, m2))
    assert m3.GetNV() == m1.GetNV() + m2.GetNV(), "merge did not work"
    print("test (3) ok")
    # m3.Print('merged_mesh2.mesh')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':
        import mfem.par as mfem
    else:
        import mfem.ser as mfem

    run_test1(mfem)
    run_test2(mfem)
    run_test3(mfem)
