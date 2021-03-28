import mfem.ser as mfem
import io

mesh = mfem.Mesh(3)
fec = mfem.H1_FECollection(1)
fes = mfem.FiniteElementSpace(mesh, fec)
x = mfem.GridFunction(fes)
x.Assign(0.0)

o = io.StringIO()
l1 = mesh.WriteToStream(o)
l2 = x.WriteToStream(o)

v = mfem.Vector([1,2,3])
l3 = v.WriteToStream(o)

print("result length: ", l1, " ", l2)
print('result: ', o.getvalue())
