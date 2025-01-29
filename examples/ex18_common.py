'''
    MFEM Example 18 - Serial/Parallel Shared Code

    This is a python translation of ex18.hpp


'''
import numpy as np
from mfem import mfem_mode

if mfem_mode is None or mfem_mode == 'serial':
    import mfem.ser as mfem
else:
    import mfem.par as mfem

from os.path import expanduser, join, dirname


class DGHyperbolicConservationLaws(mfem.TimeDependentOperator):
    def __init__(self, vfes_, formIntegrator_, preassembleWeakDivergence=True):

        super(DGHyperbolicConservationLaws, self).__init__(
            vfes_.GetTrueVSize())
        self.num_equations = formIntegrator_.num_equations
        self.vfes = vfes_
        self.dim = vfes_.GetMesh().SpaceDimension()
        self.formIntegrator = formIntegrator_

        self.z = mfem.Vector(vfes_.GetTrueVSize())

        self.weakdiv = None
        self.max_char_speed = None

        self.ComputeInvMass()

        if mfem_mode == 'serial':
            self.nonlinearForm = mfem.NonlinearForm(self.vfes)
        else:
            if isinstance(self.vfes, mfem.ParFiniteElementSpace):
                self.nonlinearForm = mfem.ParNonlinearForm(self.vfes)
            else:
                self.nonlinearForm = mfem.NonlinearForm(self.vfes)

        if preassembleWeakDivergence:
            self.ComputeWeakDivergence()
        else:
            self.nonlinearForm.AddDomainIntegrator(self.formIntegrator)

        self.nonlinearForm.AddInteriorFaceIntegrator(self.formIntegrator)
        self.nonlinearForm.UseExternalIntegrators()

    def GetMaxCharSpeed(self):
        return self.max_char_speed

    def ComputeInvMass(self):
        inv_mass = mfem.InverseIntegrator(mfem.MassIntegrator())

        self.invmass = [None]*self.vfes.GetNE()
        for i in range(self.vfes.GetNE()):
            dof = self.vfes.GetFE(i).GetDof()
            self.invmass[i] = mfem.DenseMatrix(dof)
            inv_mass.AssembleElementMatrix(self.vfes.GetFE(i),
                                           self.vfes.GetElementTransformation(
                                               i),
                                           self.invmass[i])

    def ComputeWeakDivergence(self):
        weak_div = mfem.TransposeIntegrator(mfem.GradientIntegrator())

        weakdiv_bynodes = mfem.DenseMatrix()

        self.weakdiv = [None]*self.vfes.GetNE()

        for i in range(self.vfes.GetNE()):
            dof = self.vfes.GetFE(i).GetDof()
            weakdiv_bynodes.SetSize(dof, dof*self.dim)
            weak_div.AssembleElementMatrix2(self.vfes.GetFE(i),
                                            self.vfes.GetFE(i),
                                            self.vfes.GetElementTransformation(
                                                i),
                                            weakdiv_bynodes)
            self.weakdiv[i] = mfem.DenseMatrix()
            self.weakdiv[i].SetSize(dof, dof*self.dim)

            # Reorder so that trial space is ByDim.
            # This makes applying weak divergence to flux value simpler.
            for j in range(dof):
                for d in range(self.dim):
                    self.weakdiv[i].SetCol(
                        j*self.dim + d, weakdiv_bynodes.GetColumn(d*dof + j))

    def Mult(self, x, y):
        # 0. Reset wavespeed computation before operator application.
        self.formIntegrator.ResetMaxCharSpeed()

        # 1. Apply Nonlinear form to obtain an auxiliary result
        #         z = - <F̂(u_h,n), [[v]]>_e
        #    If weak-divergence is not preassembled, we also have weak-divergence
        #         z = - <F̂(u_h,n), [[v]]>_e + (F(u_h), ∇v)
        self.nonlinearForm.Mult(x, self.z)
        #print("!!!!", self.weakdiv)
        if self.weakdiv is not None:  # if weak divergence is pre-assembled
            # Apply weak divergence to F(u_h), and inverse mass to z_loc + weakdiv_loc

            current_state = mfem.Vector()     # view of current state at a node
            current_flux = mfem.DenseMatrix()  # flux of current state

            # element flux value. Whose column is ordered by dim.
            flux = mfem.DenseMatrix()
            # view of current states in an element, dof x num_eq
            current_xmat = mfem.DenseMatrix()
            # view of element auxiliary result, dof x num_eq
            current_zmat = mfem.DenseMatrix()
            current_ymat = mfem.DenseMatrix()  # view of element result, dof x num_eq

            fluxFunction = self.formIntegrator.GetFluxFunction()

            xval = mfem.Vector()
            zval = mfem.Vector()
            flux_vec = mfem.Vector()

            for i in range(self.vfes.GetNE()):
                Tr = self.vfes.GetElementTransformation(i)
                dof = self.vfes.GetFE(i).GetDof()
                vdofs = mfem.intArray(self.vfes.GetElementVDofs(i))

                x.GetSubVector(vdofs, xval)
                current_xmat.UseExternalData(
                    xval.GetData(), dof, self.num_equations)

                #
                # Python Note:
                #    C++ code access to array data with offset is done bu GetData() + offset
                #    In Python, the same can be done by using numpy array generated from Vector::GetDataArray(),
                #
                #       array = vec.GetDataArray()
                #       new_data_pointer = mfem.Vector(array[10:]).GetData()
                #
                #    note that the above does not work if mfem.Vector is replaced by mfem.DenseMatrix
                #    This is because, while MFEM stores data in colume-major, Python numpy store raw-major.
                #

                flux.SetSize(self.num_equations, self.dim*dof)
                flux_vec = mfem.Vector(
                    flux.GetData(), self.num_equations*self.dim*dof)
                data = flux_vec.GetDataArray()

                for j in range(dof):  # compute flux for all nodes in the element
                    current_xmat.GetRow(j, current_state)

                    data_ptr = mfem.Vector(
                        data[self.num_equations*self.dim*j:]).GetData()
                    current_flux = mfem.DenseMatrix(data_ptr,
                                                    self.num_equations, dof)
                    fluxFunction.ComputeFlux(current_state, Tr, current_flux)

                # Compute weak-divergence and add it to auxiliary result, z
                # Recalling that weakdiv is reordered by dim, we can apply
                # weak-divergence to the transpose of flux.
                self.z.GetSubVector(vdofs, zval)
                current_zmat.UseExternalData(
                    zval.GetData(), dof, self.num_equations)
                mfem.AddMult_a_ABt(1.0, self.weakdiv[i], flux, current_zmat)

                # Apply inverse mass to auxiliary result to obtain the final result
                current_ymat.SetSize(dof, self.num_equations)
                mfem.Mult(self.invmass[i], current_zmat, current_ymat)
                y.SetSubVector(vdofs, current_ymat.GetData())

        else:

            # Apply block inverse mass
            zval = mfem.Vector()  # / z_loc, dof*num_eq

            # view of element auxiliary result, dof x num_eq
            current_zmat = mfem.DenseMatrix()
            current_ymat = mfem.DenseMatrix()  # view of element result, dof x num_eq

            for i in range(self.vfes.GetNE()):
                dof = self.vfes.GetFE(i).GetDof()
                vdofs = mfem.intArray(self.vfes.GetElementVDofs(i))
                self.z.GetSubVector(vdofs, zval)
                current_zmat.UseExternalData(
                    zval.GetData(), dof, self.num_equations)
                current_ymat.SetSize(dof, self.num_equations)
                mfem.Mult(self.invmass[i], current_zmat, current_ymat)
                y.SetSubVector(vdofs, current_ymat.GetData())

        self.max_char_speed = self.formIntegrator.GetMaxCharSpeed()

    def Update(self):
        self.nonlinearForm.Update()
        height = self.nonlinearForm.Height()
        width = height
        self.z.SetSize(height)

        ComputeInvMass()
        if self.weakdiv is None:
            self.ComputeWeakDivergence()


def GetMovingVortexInit(radius, Minf, beta, gas_constant, specific_heat_ratio):
    def func(x, y):
        xc = 0.0
        yc = 0.0
        # Nice units
        vel_inf = 1.
        den_inf = 1.

        pres_inf = (den_inf / specific_heat_ratio) * \
            (vel_inf / Minf) * (vel_inf / Minf)
        temp_inf = pres_inf / (den_inf * gas_constant)

        r2rad = 0.0
        r2rad += (x[0] - xc) * (x[0] - xc)
        r2rad += (x[1] - yc) * (x[1] - yc)
        r2rad /= (radius * radius)

        shrinv1 = 1.0 / (specific_heat_ratio - 1.)

        velX = vel_inf * \
            (1 - beta * (x[1] - yc) / radius * np.exp(-0.5 * r2rad))
        velY = vel_inf * beta * (x[0] - xc) / radius * np.exp(-0.5 * r2rad)
        vel2 = velX * velX + velY * velY

        specific_heat = gas_constant * specific_heat_ratio * shrinv1

        temp = temp_inf - (0.5 * (vel_inf * beta) *
                           (vel_inf * beta) / specific_heat * np.exp(-r2rad))

        den = den_inf * (temp/temp_inf)**shrinv1
        pres = den * gas_constant * temp
        energy = shrinv1 * pres / den + 0.5 * vel2

        y[0] = den
        y[1] = den * velX
        y[2] = den * velY
        y[3] = den * energy

    return func


def EulerMesh(meshfile, problem):
    if meshfile == '':
        if problem in (1, 2, 3):
            meshfile = "periodic-square.mesh"

        elif problem == 4:
            meshfile = "periodic-segment.mesh"

        else:
            assert False, "Default mesh file not given for problem = " + \
                str(problem)

        meshfile = expanduser(join(dirname(__file__), '..', 'data', meshfile))

    return mfem.Mesh(meshfile, 1, 1)

# Initial condition


def EulerInitialCondition(problem, specific_heat_ratio, gas_constant):

    if problem == 1:
        # fast moving vortex
        func = GetMovingVortexInit(0.2, 0.5, 1. / 5., gas_constant,
                                   specific_heat_ratio)
        return mfem.jit.vector(vdim=4, interface="c++")(func)

    elif problem == 2:
        # slow moving vortex
        func = GetMovingVortexInit(0.2, 0.05, 1. / 50., gas_constant,
                                   specific_heat_ratio)
        return mfem.jit.vector(vdim=4, interface="c++")(func)

    elif problem == 3:
        # moving sine wave
        @ mfem.jit.vector(vdim=4, interface="c++")
        def func(x, y):
            assert len(x) > 2, "2D is not supportd for this probl"
            density = 1.0 + 0.2 * np.sin(np.pi*(x[0]+x[1]))
            velocity_x = 0.7
            velocity_y = 0.3
            pressure = 1.0
            energy = (pressure / (1.4 - 1.0) +
                      density * 0.5 * (velocity_x * velocity_x + velocity_y * velocity_y))

            y[0] = density
            y[1] = density * velocity_x
            y[2] = density * velocity_y
            y[3] = energy

        return func

    elif problem == 4:
        @ mfem.jit.vector(vdim=3, interface="c++")
        def func(x, y):
            density = 1.0 + 0.2 * np.sin(np.pi * 2 * x[0])
            velocity_x = 1.0
            pressure = 1.0
            energy = pressure / (1.4 - 1.0) + density * \
                0.5 * (velocity_x * velocity_x)

            y[0] = density
            y[1] = density * velocity_x
            y[2] = energy

        return func

    else:
        assert False, "Problem Undefined"
