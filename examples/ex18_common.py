'''
    MFEM Example 18 - Serial/Parallel Shared Code

    This is a python translation of ex18.hpp

    note: following variabls are set from ex18 or ex18p
        problem
        num_equation
        max_char_speed
        specific_heat_ratio;
        gas_constant;

'''
import numpy as np
from mfem import mfem_mode

if mfem_mode is None or mfem_mode == 'serial':
    import mfem.ser as mfem
else:
    import mfem.par as mfem

num_equation = 0
specific_heat_ratio = 0
gas_constant = 0
problem = 0
max_char_speed = 0

class DGHyperbolicConservationLaws(mfem.TimeDependentOperator):
    def __init__(self, vfes_, formIntegrator_, preassembleWeakDivergence=True):
        
        super(DGHyperbolicConservationLaws, self).__init__(vfes_.GetTrueVsize())
        self.num_equations = formIntegrator_.num_equations
        self.vfes = vfes_
        self.dim = vfes_.GetMesh().SpaceDimension()
        self.formIntegrator = formIntegrator_

        self.z = mfem.Vector(vfes_.GetTrueVSize())
        self.weakdiv = None

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
        self.nonlinearForm.UseExternalIntegrator()

    def ComputeInvMass(self):
        inv_mass = mfem.InverseIntegrator(mfem.MassIntegrator())

        self.invmass = [None]*self.vfes.GetNE()
        for i in range(self.vfes.GetNE()):
            dof = vfes.GetFE(i).GetDof()
            invmass[i] = mfem.DenseMatrix(dof)
            inv_mass.AssembleElementMatrix(self.vfes.GetFE(i),
                                           self.vfes.GetElementTransformation(i),
                                           invmass[i])
        

    def ComputeWeakDivergence(self):
        weak_div = mfem.TransposeIntegrator(mfem.GradientIntegrator())
        
        weakdiv_bynodes = mfem.DenseMatrix()

        self.weakdiv = [None]*self.vfes.GetNE()

        for i in range(self.vfes.GetNE()):
            dof = vfes.GetFE(i).GetDof()
            weakdiv_bynodes.SetSize(dof, dof*self.dim)
            weak_div.AssembleElementMatrix2(self.vfes.GetFE(i),
                                            self.vfes.GetFE(i),
                                            self.vfes.GetElementTransformation(i),
                                            weakdiv_bynodes);
            self.weakdiv[i] = mfem.DenseMatrix()
            self.weakdiv[i].SetSize(dof, dof*dim)
            
            # Reorder so that trial space is ByDim.
            # This makes applying weak divergence to flux value simpler.
            for j in range(dof):
                for d in range(self.dim):
                    self.weakdiv[i].SetCol(j*dim + d, weakdiv_bynodes.GetColumn(d*dof + j))


    def Mult(x, y):
        # 0. Reset wavespeed computation before operator application.
        self.formIntegrator.ResetMaxCharSpeed()

        
        # 1. Apply Nonlinear form to obtain an auxiliary result
        #         z = - <F̂(u_h,n), [[v]]>_e
        #    If weak-divergence is not preassembled, we also have weak-divergence
        #         z = - <F̂(u_h,n), [[v]]>_e + (F(u_h), ∇v)
        self.nonlinearForm.>Mult(x, z)
        
        if self.weakdiv is not None: # if weak divergence is pre-assembled
             # Apply weak divergence to F(u_h), and inverse mass to z_loc + weakdiv_loc
            
             current_state = mfem.Vector()     # view of current state at a node
             current_flux = mfem.DenseMatrix() # flux of current state

             flux = mfem.DenseMatrix()         # element flux value. Whose column is ordered by dim.
             current_xmat = mfem.DenseMatrix() # view of current states in an element, dof x num_eq
             current_zmat = mfem.DenseMatrix() # view of element auxiliary result, dof x num_eq
             current_ymat = mfem.DenseMatrix() # view of element result, dof x num_eq

             fluxFunction = formIntegrator.GetFluxFunction()
             
             vdofs = mfem.intArray()
             xval = mfem.Vector()
             zval = mfem.Vector()
             
             for i in range(self.vfes.GetNE()):
                 Tr = self.vfes.GetElementTransformation(i)
                 dof = self.vfes.GetFE(i).GetDof()
                 self.vfes.GetElementVDofs(i, vdofs)
         x.GetSubVector(vdofs, xval);
         current_xmat.UseExternalData(xval.GetData(), dof, num_equations);
         flux.SetSize(num_equations, dim*dof);
         for (int j=0; j<dof; j++) // compute flux for all nodes in the element
         {
            current_xmat.GetRow(j, current_state);
            current_flux.UseExternalData(flux.GetData() + num_equations*dim*j,
                                         num_equations, dof);
            fluxFunction.ComputeFlux(current_state, *Tr, current_flux);
         }
         // Compute weak-divergence and add it to auxiliary result, z
         // Recalling that weakdiv is reordered by dim, we can apply
         // weak-divergence to the transpose of flux.
         z.GetSubVector(vdofs, zval);
         current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
         mfem::AddMult_a_ABt(1.0, weakdiv[i], flux, current_zmat);
         // Apply inverse mass to auxiliary result to obtain the final result
         current_ymat.SetSize(dof, num_equations);
         mfem::Mult(invmass[i], current_zmat, current_ymat);
         y.SetSubVector(vdofs, current_ymat.GetData());
      }
   }
   else
   {
      // Apply block inverse mass
      Vector zval; // z_loc, dof*num_eq

      DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
      DenseMatrix current_ymat; // view of element result, dof x num_eq
      Array<int> vdofs;
      for (int i=0; i<vfes.GetNE(); i++)
      {
         int dof = vfes.GetFE(i)->GetDof();
         vfes.GetElementVDofs(i, vdofs);
         z.GetSubVector(vdofs, zval);
         current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
         current_ymat.SetSize(dof, num_equations);
         mfem::Mult(invmass[i], current_zmat, current_ymat);
         y.SetSubVector(vdofs, current_ymat.GetData());
      }
   }
   max_char_speed = formIntegrator->GetMaxCharSpeed();
}
                    


class FE_Evolution(mfem.TimeDependentOperator):
    def __init__(self, vfes, A, A_flux):
        self.dim = vfes.GetFE(0).GetDim()
        self.vfes = vfes
        self.A = A
        self.Aflux = A_flux
        self.Me_inv = mfem.DenseTensor(vfes.GetFE(0).GetDof(),
                                       vfes.GetFE(0).GetDof(),
                                       vfes.GetNE())

        self.state = mfem.Vector(num_equation)
        self.f = mfem.DenseMatrix(num_equation, self.dim)
        self.flux = mfem.DenseTensor(vfes.GetNDofs(), self.dim, num_equation)


        dof = vfes.GetFE(0).GetDof()
        Me = mfem.DenseMatrix(dof)
        inv = mfem.DenseMatrixInverse(Me)
        mi = mfem.MassIntegrator()
        for i in range(vfes.GetNE()):
            mi.AssembleElementMatrix(vfes.GetFE(
                i), vfes.GetElementTransformation(i), Me)
            inv.Factor()
            inv.GetInverseMatrix(self.Me_inv(i))
        super(FE_Evolution, self).__init__(A.Height())

    def GetFlux(self, x, flux):
        state = self.state
        dof = self.flux.SizeI()
        dim = self.flux.SizeJ()

        flux_data = []
        for i in range(dof):
            for k in range(num_equation):
                self.state[k] = x[i, k]
            ComputeFlux(state, dim, self.f)

            flux_data.append(self.f.GetDataArray().transpose().copy())
            # flux[i].Print()
            # print(self.f.GetDataArray())
            # for d in range(dim):
            #    for k in range(num_equation):
            #        flux[i, d, k] = self.f[k, d]

            mcs = ComputeMaxCharSpeed(state, dim)
            if (mcs > globals()['max_char_speed']):
                globals()['max_char_speed'] = mcs

        flux.Assign(np.stack(flux_data))
        # print("max char speed", globals()['max_char_speed'])

    def Mult(self, x, y):
        globals()['max_char_speed'] = 0.
        num_equation = globals()['num_equation']
        # 1. Create the vector z with the face terms -<F.n(u), [w]>.
        self.A.Mult(x, self.z)

        # 2. Add the element terms.
        # i.  computing the flux approximately as a grid function by interpolating
        #     at the solution nodes.
        # ii. multiplying this grid function by a (constant) mixed bilinear form for
        #     each of the num_equation, computing (F(u), grad(w)) for each equation.

        xmat = mfem.DenseMatrix(
            x.GetData(), self.vfes.GetNDofs(), num_equation)
        self.GetFlux(xmat, self.flux)

        for k in range(num_equation):
            fk = mfem.Vector(self.flux[k].GetData(),
                             self.dim * self.vfes.GetNDofs())
            o = k * self.vfes.GetNDofs()
            zk = self.z[o: o+self.vfes.GetNDofs()]
            self.Aflux.AddMult(fk, zk)

        # 3. Multiply element-wise by the inverse mass matrices.
        zval = mfem.Vector()
        vdofs = mfem.intArray()
        dof = self.vfes.GetFE(0).GetDof()
        zmat = mfem.DenseMatrix()
        ymat = mfem.DenseMatrix(dof, num_equation)

        for i in range(self.vfes.GetNE()):
            # Return the vdofs ordered byNODES
            vdofs = mfem.intArray(self.vfes.GetElementVDofs(i))
            self.z.GetSubVector(vdofs, zval)
            zmat.UseExternalData(zval.GetData(), dof, num_equation)
            mfem.Mult(self.Me_inv[i], zmat, ymat)
            y.SetSubVector(vdofs, ymat.GetData())


class DomainIntegrator(mfem.PyBilinearFormIntegrator):
    def __init__(self, dim):
        num_equation = globals()['num_equation']
        self.flux = mfem.DenseMatrix(num_equation, dim)
        self.shape = mfem.Vector()
        self.dshapedr = mfem.DenseMatrix()
        self.dshapedx = mfem.DenseMatrix()
        super(DomainIntegrator, self).__init__()

    def AssembleElementMatrix2(self, trial_fe, test_fe, Tr, elmat):
        # Assemble the form (vec(v), grad(w))

        # Trial space = vector L2 space (mesh dim)
        # Test space  = scalar L2 space

        dof_trial = trial_fe.GetDof()
        dof_test = test_fe.GetDof()
        dim = trial_fe.GetDim()

        self.shape.SetSize(dof_trial)
        self.dshapedr.SetSize(dof_test, dim)
        self.dshapedx.SetSize(dof_test, dim)

        elmat.SetSize(dof_test, dof_trial * dim)
        elmat.Assign(0.0)

        maxorder = max(trial_fe.GetOrder(), test_fe.GetOrder())
        intorder = 2 * maxorder
        ir = mfem.IntRules.Get(trial_fe.GetGeomType(), intorder)

        for i in range(ir.GetNPoints()):
            ip = ir.IntPoint(i)

            # Calculate the shape functions
            trial_fe.CalcShape(ip, self.shape)
            self.shape *= ip.weight

            # Compute the physical gradients of the test functions
            Tr.SetIntPoint(ip)
            test_fe.CalcDShape(ip, self.dshapedr)
            mfem.Mult(self.dshapedr, Tr.AdjugateJacobian(), self.dshapedx)

            for d in range(dim):
                for j in range(dof_test):
                    for k in range(dof_trial):
                        elmat[j, k + d * dof_trial] += self.shape[k] * \
                            self.dshapedx[j, d]


class FaceIntegrator(mfem.PyNonlinearFormIntegrator):
    def __init__(self, rsolver, dim):
        self.rsolver = rsolver
        self.shape1 = mfem.Vector()
        self.shape2 = mfem.Vector()
        self.funval1 = mfem.Vector(num_equation)
        self.funval2 = mfem.Vector(num_equation)
        self.nor = mfem.Vector(dim)
        self.fluxN = mfem.Vector(num_equation)
        self.eip1 = mfem.IntegrationPoint()
        self.eip2 = mfem.IntegrationPoint()
        super(FaceIntegrator, self).__init__()

        self.fluxNA = np.atleast_2d(self.fluxN.GetDataArray())

    def AssembleFaceVector(self, el1, el2, Tr, elfun, elvect):
        num_equation = globals()['num_equation']
        # Compute the term <F.n(u),[w]> on the interior faces.
        dof1 = el1.GetDof()
        dof2 = el2.GetDof()

        self.shape1.SetSize(dof1)
        self.shape2.SetSize(dof2)

        elvect.SetSize((dof1 + dof2) * num_equation)
        elvect.Assign(0.0)

        elfun1_mat = mfem.DenseMatrix(elfun.GetData(), dof1, num_equation)
        elfun2_mat = mfem.DenseMatrix(
            elfun[dof1*num_equation:].GetData(), dof2, num_equation)

        elvect1_mat = mfem.DenseMatrix(elvect.GetData(), dof1, num_equation)
        elvect2_mat = mfem.DenseMatrix(
            elvect[dof1*num_equation:].GetData(), dof2, num_equation)

        # Integration order calculation from DGTraceIntegrator
        if (Tr.Elem2No >= 0):
            intorder = (min(Tr.Elem1.OrderW(), Tr.Elem2.OrderW()) +
                        2*max(el1.GetOrder(), el2.GetOrder()))
        else:
            intorder = Tr.Elem1.OrderW() + 2*el1.GetOrder()

        if (el1.Space() == mfem.FunctionSpace().Pk):
            intorder += 1

        ir = mfem.IntRules.Get(Tr.GetGeometryType(), int(intorder))

        mat1A = elvect1_mat.GetDataArray()
        mat2A = elvect2_mat.GetDataArray()
        shape1A = np.atleast_2d(self.shape1.GetDataArray())
        shape2A = np.atleast_2d(self.shape2.GetDataArray())

        for i in range(ir.GetNPoints()):
            ip = ir.IntPoint(i)
            Tr.Loc1.Transform(ip, self.eip1)
            Tr.Loc2.Transform(ip, self.eip2)

            # Calculate basis functions on both elements at the face
            el1.CalcShape(self.eip1, self.shape1)
            el2.CalcShape(self.eip2, self.shape2)

            # Interpolate elfun at the point
            elfun1_mat.MultTranspose(self.shape1, self.funval1)
            elfun2_mat.MultTranspose(self.shape2, self.funval2)
            Tr.Face.SetIntPoint(ip)

            # Get the normal vector and the flux on the face

            mfem.CalcOrtho(Tr.Face.Jacobian(), self.nor)

            mcs = self.rsolver.Eval(
                self.funval1, self.funval2, self.nor, self.fluxN)

            # Update max char speed
            if mcs > globals()['max_char_speed']:
                globals()['max_char_speed'] = mcs

            self.fluxN *= ip.weight

            #
            mat1A -= shape1A.transpose().dot(self.fluxNA)
            mat2A += shape2A.transpose().dot(self.fluxNA)
            '''
            for k in range(num_equation):
                for s in range(dof1):
                    elvect1_mat[s, k] -= self.fluxN[k] * self.shape1[s]
                for s in range(dof2):
                    elvect2_mat[s, k] += self.fluxN[k] * self.shape2[s]
            '''


class RiemannSolver(object):
    def __init__(self):
        num_equation = globals()['num_equation']
        self.flux1 = mfem.Vector(num_equation)
        self.flux2 = mfem.Vector(num_equation)

    def Eval(self, state1, state2, nor, flux):

        # NOTE: nor in general is not a unit normal
        dim = nor.Size()

        assert StateIsPhysical(state1, dim), ""
        assert StateIsPhysical(state2, dim), ""

        maxE1 = ComputeMaxCharSpeed(state1, dim)
        maxE2 = ComputeMaxCharSpeed(state2, dim)
        maxE = max(maxE1, maxE2)

        ComputeFluxDotN(state1, nor, self.flux1)
        ComputeFluxDotN(state2, nor, self.flux2)

        # normag = np.sqrt(np.sum(nor.GetDataArray()**2))
        normag = nor.Norml2()

        '''
        for i in range(num_equation):
            flux[i] = (0.5 * (self.flux1[i] + self.flux2[i])
                       - 0.5 * maxE * (state2[i] - state1[i]) * normag)
        '''
        f = (0.5 * (self.flux1.GetDataArray() + self.flux2.GetDataArray())
             - 0.5 * maxE * (state2.GetDataArray() - state1.GetDataArray()) * normag)
        flux.Assign(f)

        return maxE


def StateIsPhysical(state, dim):
    specific_heat_ratio = globals()["specific_heat_ratio"]

    den = state[0]
    # den_vel = state.GetDataArray()[1:1+dim]
    den_energy = state[1 + dim]

    if (den < 0):
        print("Negative density: " + str(state.GetDataArray()))
        return False
    if (den_energy <= 0):
        print("Negative energy: " + str(state.GetDataArray()))
        return False

    # den_vel2 = np.sum(den_vel**2)/den
    den_vel2 = (state[1:1+dim].Norml2())**2/den
    pres = (specific_heat_ratio - 1.0) * (den_energy - 0.5 * den_vel2)
    if (pres <= 0):
        print("Negative pressure: " + str(state.GetDataArray()))
        return False
    return True


class InitialCondition(mfem.VectorPyCoefficient):
    def __init__(self, dim):
        mfem.VectorPyCoefficient.__init__(self, dim)

    def EvalValue(self, x):
        dim = x.shape[0]
        assert dim == 2, ""
        problem = globals()['problem']
        if (problem == 1):
            # "Fast vortex"
            radius = 0.2
            Minf = 0.5
            beta = 1. / 5.
        elif (problem == 2):
            # "Slow vortex"
            radius = 0.2
            Minf = 0.05
            beta = 1. / 50.
        else:
            assert False, "Cannot recognize problem. Options are: 1 - fast vortex, 2 - slow vortex"

        xc = 0.0
        yc = 0.0
        # Nice units
        vel_inf = 1.
        den_inf = 1.

        specific_heat_ratio = globals()["specific_heat_ratio"]
        gas_constant = globals()["gas_constant"]

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

        y = np.array([den, den * velX, den * velY, den * energy])
        return y


def ComputePressure(state, dim):
    den = state[0]
    # den_vel = state.GetDataArray()[1:1+dim]
    den_energy = state[1 + dim]

    specific_heat_ratio = globals()["specific_heat_ratio"]
    # den_vel2 = np.sum(den_vel**2)/den
    den_vel2 = (state[1:1+dim].Norml2())**2/den
    pres = (specific_heat_ratio - 1.0) * (den_energy - 0.5 * den_vel2)

    return pres


def ComputeFlux(state, dim, flux):
    den = state[0]
    den_vel = state.GetDataArray()[1:1+dim]
    den_energy = state[1 + dim]

    assert StateIsPhysical(state, dim), ""

    pres = ComputePressure(state, dim)

    den_vel2 = np.atleast_2d(den_vel)
    fluxA = flux.GetDataArray()
    fluxA[0, :] = den_vel
    fluxA[1:1+dim, :] = den_vel2.transpose().dot(den_vel2) / den
    for d in range(dim):
        fluxA[1+d, d] += pres

    H = (den_energy + pres) / den
    flux.GetDataArray()[1+dim, :] = den_vel * H


def ComputeFluxDotN(state, nor, fluxN):
    # NOTE: nor in general is not a unit normal
    dim = nor.Size()
    nor = nor.GetDataArray()
    fluxN = fluxN.GetDataArray()

    den = state[0]
    den_vel = state.GetDataArray()[1:1+dim]
    den_energy = state[1 + dim]

    assert StateIsPhysical(state, dim), ""

    pres = ComputePressure(state, dim)

    den_velN = den_vel.dot(nor)

    fluxN[0] = den_velN
    fluxN[1:1+dim] = den_velN * den_vel / den + pres * nor

    H = (den_energy + pres) / den
    fluxN[1+dim] = den_velN * H


def ComputeMaxCharSpeed(state, dim):
    specific_heat_ratio = globals()["specific_heat_ratio"]

    den = state[0]
    den_vel2 = (state[1:1+dim].Norml2())**2/den
    pres = ComputePressure(state, dim)

    sound = np.sqrt(specific_heat_ratio * pres / den)
    vel = np.sqrt(den_vel2 / den)

    return vel + sound


std: : function < void(const Vector&, Vector&) > GetMovingVortexInit(
   const real_t radius, const real_t Minf, const real_t beta,
   const real_t gas_constant, const real_t specific_heat_ratio)
{
   return [specific_heat_ratio,
           gas_constant, Minf, radius, beta](const Vector & x, Vector & y)
   {
      MFEM_ASSERT(x.Size() == 2, "");

      const real_t xc = 0.0, yc = 0.0;

      // Nice units
      const real_t vel_inf = 1.;
      const real_t den_inf = 1.;

      // Derive remainder of background state from this and Minf
      const real_t pres_inf = (den_inf / specific_heat_ratio) *
                              (vel_inf / Minf) * (vel_inf / Minf);
      const real_t temp_inf = pres_inf / (den_inf * gas_constant);

def GetMovingVertexInit(radius, Minf, beta, gas_constant, specific_heat_ratio):
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

def EulerMesh(meshfile, problem)
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
