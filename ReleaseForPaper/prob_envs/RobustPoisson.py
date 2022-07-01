from gym import spaces
import numpy as np
from prob_envs.Poisson import Poisson
from utils.ReentrantCorner import *
import os
from os.path import expanduser, join
from mfem.ser import intArray

class RobustAnglePoisson(Poisson):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)  
        problem_type = kwargs.get('problem_type','lshaped')
        if (problem_type == 'lshaped'):
            delattr(self, 'BC')
            self.BC = mfem.NumbaFunction(ReentrantCornerExact, 2, True).GenerateCoefficient()
        elif (problem_type == 'wavefront'):
            print("*** Should not be using RobustAnglePoisson class with wavefront problem - exiting ***")
            exit()
        else:
            delattr(self, 'BC')
            self.BC = mfem.ConstantCoefficient(0.0)
        self.num_unif_ref = kwargs.get('num_unif_ref',1)        
        self.Lshapedmeshfile = expanduser(join(os.path.dirname(__file__), '../..', 'data', 'l-shape-benchmark.mesh'))
        self.circlemeshfile = expanduser(join(os.path.dirname(__file__), '../..', 'data', 'circle_3_4.mesh'))

        self.angle = kwargs.get('angle_lower', np.pi * 0.5)
        self.angle_lower = kwargs.get('angle_lower', np.pi * 0.25)
        self.angle_upper = kwargs.get('angle_upper', np.pi * 0.75)
        if (problem_type == 'lshaped'):
            self.set_angle(self.angle)

    def set_angle(self, angle):
        # print("Setting env angle to ", self.angle, flush=True)
        self.angle = angle
        self.BC.SetTime(self.angle)

        if self.mesh_name =='l-shape-benchmark.mesh':
            self.initial_mesh = ReentrantCornerMesh(self.angle, self.Lshapedmeshfile)
        elif self.mesh_name =='circle_3_4.mesh':
            self.initial_mesh = ReentrantCornerMesh(self.angle, self.circlemeshfile)
        else:
            print("Something went wrong - see set_angle routine; exiting")
            exit()
        for _ in range(self.num_unif_ref):
            self.initial_mesh.UniformRefinement()
        self.initial_mesh.EnsureNCMesh()

    def reset(self, random_angle=True):
        if random_angle:
            angle = np.random.uniform(self.angle_lower, self.angle_upper, 1).item()
            self.set_angle(angle)
        return super().reset()


class hpPoisson(RobustAnglePoisson):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.action_space = spaces.Box(low=0.0, high=0.999, shape=(2,), dtype=np.float32)

    # p-refine for 'max' refinement strategy: 
    #   inputs theta and rho are in [0,1]
    #   1) collect elements with errors >= theta * (max error)
    #   2) from that collection, p-refine if and only if error >= rho * (max error)
    #   effectively: 
    #      if (theta <= rho): p-refine all elements with error >= rho * (max error)
    #      else:              p-refine all elements with error >= theta * (max error)  
    def Prefine(self, theta, rho):   
        mark_to_p_refine = []
        threshold = theta * np.max(self.errors)
        for i in range(self.mesh.GetNE()):
            if threshold >= self.errors[i]:
                mark_to_p_refine.append((i, self.errors[i]))
        mark_to_p_refine.sort(key=lambda x:x[1], reverse=True)
        for i in range(len(mark_to_p_refine)):
            if mark_to_p_refine[i][1] > rho * np.max(self.errors):
                current_element = mark_to_p_refine[i][0]
                current_order = self.fespace.GetElementOrder(current_element)
                self.fespace.SetElementOrder(current_element, current_order + 1)
        
        self.fespace.Update(False)
        self.x.Update()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        self.a.Update()
        self.b.Update()

    def CloseMesh(self, delta_p = 1):
        # Loop through all elements in mesh until the maximum difference in polynomial
        # orders across all edges is no more than delta_p
        neighbor_table = self.mesh.ElementToElementTable()
        while True:
            mesh_closed = True
            elements_to_p_refine = []
            for i in range(self.mesh.GetNE()):
                neighbor_row = neighbor_table.GetRow(i)
                row_size = neighbor_table.RowSize(i)
                neighbor_array = intArray(row_size)
                neighbor_array.Assign(neighbor_row)
                for l in range(row_size):
                    neighbor_order = self.fespace.GetElementOrder(neighbor_array[l])
                    if neighbor_order - self.fespace.GetElementOrder(i) > delta_p:
                        elements_to_p_refine.append(i)
                        mesh_closed = False
            p_refine_elements = np.unique(elements_to_p_refine).tolist()
            for k in range(len(p_refine_elements)):
                current_element = p_refine_elements[k]
                current_order = self.fespace.GetElementOrder(current_element)
                self.fespace.SetElementOrder(current_element, current_order + 1)
            if mesh_closed:
                break

        self.fespace.Update(False)
        self.x.Update()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        # self.fespace.UpdatesFinished()
        self.a.Update()
        self.b.Update()

    # overriding UpdateMesh in Poisson (grand)-parent class
    def UpdateMesh(self, action):
        action = np.clip(action, 0.0, 1.0)
        theta = action[0].item() 
        rho = action[1] * theta 
        self.Prefine(theta, rho)
        self.Refine(theta)
        self.CloseMesh()
    
    def GetElementVertices(self, k):
        Tr = self.mesh.GetElementTransformation(k)
        physical_pts = np.zeros((4,2))
        reference_pt = mfem.IntegrationPoint()
        for i in range(2):
            for j in range(2):
                reference_pt.Set(float(i),float(j),0.0,0.0)
                physical_pts[i+2*j,:] = Tr.Transform(reference_pt)
        return physical_pts


    def RenderHPmesh(self, gfname=None):
        ordersfec = mfem.L2_FECollection(0, self.dim)
        ordersfes = mfem.FiniteElementSpace(self.mesh, ordersfec)
        orders = mfem.GridFunction(ordersfes)
        for i in range(0, self.mesh.GetNE()):
            elem_dofs = 0
            elem_dofs = ordersfes.GetElementDofs(i)
            # orders[elem_dofs[0]] = self.errors[i]
            orders[elem_dofs[0]] = self.fespace.GetElementOrder(i)
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock.send_solution(self.mesh, orders)
        title = "step " + str(self.k)
        sol_sock.send_text('keys ARjlmp*******' + " window_title '" + title)
        sol_sock.send_text("valuerange 1.0 8.0 \n")
        if gfname:
            orders.Save(str(gfname), 8)  # second input is "precision"


    # action here is a single value between 0 and 1
    def expertStep(self, action):

        self.mesh.EnsureNodes()
        elements_to_h_refine = []
        elements_to_p_refine = []
        element_error_list = []

        # # here is where data is saved for the det-with-flag case.  
        # self.rows.append([action, angle, self.k, self.mesh.GetNE(), self.fespace.GetTrueVSize(), self.sum_of_dofs, self.global_error, episode_cost])

        for i in range(self.mesh.GetNE()):
            element_error_list.append((i, self.errors[i]))
        element_error_list.sort(key=lambda x:x[1], reverse=True)

        for i in range(self.mesh.GetNE()):
            vertex_coords = self.GetElementVertices(i)
            element_touching_corner = False
            curr_error = self.errors[i]
            threshold = action * np.max(self.errors)
            
            if threshold < curr_error:
                for j in range(4):
                    x = vertex_coords[j,0]
                    y = vertex_coords[j,1]
                    if abs(x) < 1e-10 and abs(y) < 1e-10:
                        element_touching_corner = True
                if(element_touching_corner):
                    elements_to_h_refine.append(i)
                else:
                    elements_to_p_refine.append(i)

        p_refine_elements = np.unique(elements_to_p_refine).tolist()
        for k in range(len(p_refine_elements)):
            current_element = p_refine_elements[k]
            current_order = self.fespace.GetElementOrder(current_element)
            self.fespace.SetElementOrder(current_element, current_order + 1)
        
        self.fespace.Update(False)
        self.x.Update()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        self.a.Update()
        self.b.Update()

        elements_to_h_refine = intArray(elements_to_h_refine)
        self.mesh.GeneralRefinement(elements_to_h_refine)
        
        self.fespace.Update(False)
        self.x.Update()
        self.x.Assign(0.0)
        self.x.ProjectBdrCoefficient(self.BC, self.ess_bdr)
        self.a.Update()
        self.b.Update()

        self.CloseMesh()

        self.AssembleAndSolve()
        self.errors = self.GetLocalErrors()
        global_error = self.GetGlobalError()
        if self.k == 1:
            cost = np.log2(global_error)
        else:
            cost = np.log2(global_error/self.global_error)
        self.global_error = global_error

        num_dofs = self.fespace.GetTrueVSize()
        self.sum_of_dofs += num_dofs

        done = False
        
        if self.sum_of_dofs > self.dof_threshold:
            cost = 0.0
            done = True

        obs = self.GetObservation()
        info = {'global_error':self.global_error, 'num_dofs':num_dofs, 'max_local_errors':np.amax(self.errors)}

        self.k += 1

        return obs, -cost, done, info