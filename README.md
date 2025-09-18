[![badge](examples/jupyter/.assets/ex1.svg)](https://colab.research.google.com/github/mfem/pymfem/blob/master/examples/jupyter/ex1.ipynb)
[![badge](examples/jupyter/.assets/ex9.svg)](https://colab.research.google.com/github/mfem/pymfem/blob/master/examples/jupyter/ex9.ipynb)

#  MFEM + PyMFEM (FEM library)

This repository provides Python binding for MFEM. MFEM is a high performance parallel finite element method (FEM) library (http://mfem.org/).

Installer (setup.py) builds both MFEM and binding together.
By default, "pip install mfem" downloads and builds the serial version of MFEM and PyMFEM.
Additionally, the installer supports building MFEM with specific options together with other external libraries, including MPI version.

## Install
### Using pip (Serial MFEM)
```shell
pip install mfem            # binary install is available only on linux platforms (Py38-312)

```

### Build with additional features (MPI, GPU, GPU-Hypre, GSLIB, SuiteSparse, libCEED, LAPACK)

The setup script accept various options. Download the package manually and run the script. Examples below downloads and build parallel version of MFEM library (linked with Metis and Hypre) and installs under <prefix>/mfem. See [INSTALL](INSTALL.md) for various other options

#### Download source and build
```shell
$ git clone https://github.com/mfem/PyMFEM.git

# Build it from local source with MPI
$ pip install ./ -C"with-parallel=Yes" --verbose
```

#### Cleaning
```shell
$ python setup.py clean --all # clean external dependencies + wrapper code
```
#### Run test
```shell
cd test
python test_examples.py -serial
```

## Usage
This example (modified from `ex1.cpp`) solves the Poisson equation,
$$\nabla \cdot (\alpha \nabla u) = f$$
in a square and plots the result using matplotlib.
Use the badge above to open this in Colab.

```python
import mfem.ser as mfem

# Create a square mesh
mesh = mfem.Mesh(10, 10, "TRIANGLE")

# Define the finite element function space
fec = mfem.H1_FECollection(1, mesh.Dimension())   # H1 order=1
fespace = mfem.FiniteElementSpace(mesh, fec)

# Define the essential dofs
ess_tdof_list = mfem.intArray()
ess_bdr = mfem.intArray([1]*mesh.bdr_attributes.Size())
fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)

# Define constants for alpha (diffusion coefficient) and f (RHS)
alpha = mfem.ConstantCoefficient(1.0)
rhs = mfem.ConstantCoefficient(1.0)

"""
Note
-----
In order to represent a variable diffusion coefficient, you
must use a numba-JIT compiled function. For example:

>>> @mfem.jit.scalar
>>> def alpha(x):
>>>     return x+1.0
"""

# Define the bilinear and linear operators
a = mfem.BilinearForm(fespace)
a.AddDomainIntegrator(mfem.DiffusionIntegrator(alpha))
a.Assemble()
b = mfem.LinearForm(fespace)
b.AddDomainIntegrator(mfem.DomainLFIntegrator(rhs))
b.Assemble()

# Initialize a gridfunction to store the solution vector
x = mfem.GridFunction(fespace)
x.Assign(0.0)

# Form the linear system of equations (AX=B)
A = mfem.OperatorPtr()
B = mfem.Vector()
X = mfem.Vector()
a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
print("Size of linear system: " + str(A.Height()))

# Solve the linear system using PCG and store the solution in x
AA = mfem.OperatorHandle2SparseMatrix(A)
M = mfem.GSSmoother(AA)
mfem.PCG(AA, M, B, X, 1, 200, 1e-12, 0.0)
a.RecoverFEMSolution(X, b, x)

# Extract vertices and solution as numpy arrays
verts = mesh.GetVertexArray()
sol = x.GetDataArray()

# Plot the solution using matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri

triang = tri.Triangulation(verts[:,0], verts[:,1])

fig, ax = plt.subplots()
ax.set_aspect('equal')
tpc = ax.tripcolor(triang, sol, shading='gouraud')
fig.colorbar(tpc)
plt.show()
```
![](https://raw.githubusercontent.com/mfem/PyMFEM/master/docs/example_image.png)


## License
PyMFEM is licensed under BSD-3 license. All new contributions
must be made under this license. See [License](LICENSE) for details.

Please refer the developers' web sites for the external libraries
* MFEM: https://mfem.org/
* Hypre: https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods
* METIS: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
* libceed: https://github.com/CEED/libCEED
* gslib: https://github.com/Nek5000/gslib
