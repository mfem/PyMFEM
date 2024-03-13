[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mfem/PyMFEM/HEAD?labpath=examples%2Fjupyter)
[![badge](examples/jupyter/ex1.svg)](https://mybinder.org/v2/gh/mfem/PyMFEM/HEAD?labpath=examples%2Fjupyter%2Fex1.ipynb)
[![badge](examples/jupyter/ex9.svg)](https://mybinder.org/v2/gh/mfem/PyMFEM/HEAD?labpath=examples%2Fjupyter%2Fex9.ipynb)



#  MFEM + PyMFEM (FEM library)

This repository repository modifies the pyMFEM setup to include the Navier miniapp. See the pyMFEM install page for full details

## Install

Run the command
`git clone https://github.com/mfem/PyMFEM.git`

Then setup a virtual environment and install the packages 
`pip install -r PyMFEM/requirements.txt`

Install PyMFEM 
`python setup.py install --with-parallel`

PyMFEM should install with `navier_solver` as a module in python that is called by `mfem.navier_solver.<function_we_want>`.

Currently only `navier_mms.cpp` is converted from MFEM. We have also included `navier_2dfocs` which simulates flow over a 2D cylinder (based on `navier_3dfocs`). Converting and other input files should be similar.


## License
PyMFEM is licensed under BSD-3.
Please refer the developers' web sites for the external libraries
* MFEM: https://mfem.org/
* Hypre: https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods
* METIS: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
* libceed: https://github.com/CEED/libCEED
* gslib: https://github.com/Nek5000/gslib
