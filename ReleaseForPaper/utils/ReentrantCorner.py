from math import atan2, sqrt, sin, cos
from numba import cfunc
import mfem.ser as mfem
import numpy as np

def ReentrantCornerMesh(omega, meshfile):
    mesh = mfem.Mesh(meshfile)
    mesh.EnsureNodes()

    order = 1
    dim = mesh.Dimension()
    fec = mfem.H1_FECollection(order, dim)
    fespace = mfem.FiniteElementSpace(mesh, fec)

    nodes = mesh.GetNodes()
    num_nodes = int(nodes.Size()/2)
    for i in range(num_nodes):
        x = nodes[2*i]
        y = nodes[2*i+1]
        theta = atan2(y, x)
        if x > 0 and abs(y) < 1e-6:
            theta = 0.0
        elif y < 0:
            theta += 2*np.pi
        delta_theta = theta/(3*np.pi/2) * (np.pi/2 - omega)
        x_tmp = x
        y_tmp = y
        x = x_tmp*cos(delta_theta) - y_tmp*sin(delta_theta)
        y = x_tmp*sin(delta_theta) + y_tmp*cos(delta_theta)
        nodes[2*i] = x
        nodes[2*i+1] = y

    return mesh

@cfunc(mfem.scalar_sig_t)
def ReentrantCornerExact(pt, omega, sdim):
    comp_omega = 2*np.pi - omega
    x = pt[0]
    y = pt[1]
    r = sqrt(x*x + y*y)
    alpha = np.pi / comp_omega
    theta = atan2(y, x)
    if x > 0 and abs(y) < 1e-6:
        theta = 0.0
    elif y < 0:
        theta += 2*np.pi
    return r**alpha * sin(alpha * theta)