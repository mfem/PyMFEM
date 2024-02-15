from mfem.par import intArray, doubleArray
import mfem.par as mfem
import os
import sys
from os.path import expanduser, join
from numpy import sin, cos, exp, sqrt, zeros, abs, pi
import numpy as np
from mpi4py import MPI


def CopyDBFIntegrators(src, dst):
    bffis = src.GetDBFI

    for i in range(bffis.Size()):
      dst.AddDomainIntegrator(bffis[i])
      
    return src, dst

def NavierSolver(mesh, order, kin_vis):
    pmesh = mfem.ParMesh(MPI.COMM_WORLD, mesh)
    mfem.IntegrationRules(0, mfem.Quadrature1D.GaussLobatto)
    vfec = mfem.H1_FECollection(order, pmesh.Dimension())
    pfec = mfem.H1_FECollection(order)
    vfes = mfem.ParFiniteElementSpace(pmesh, vfec, pmesh.Dimension())
    pfes = mfem.ParFiniteElementSpace(pmesh, pfec)

    ## Check if fully periodic mesh
    if (pmesh.bdr_attributes.Size() != 0):
        
        vel_ess_attr.SetSize(pmesh.bdr_attributes.Max())
        vel_ess_attr = 0

        pres_ess_attr.SetSize(pmesh.bdr_attributes.Max())
        pres_ess_attr = 0
    
    vfes_truevsize = vfes.GetTrueVSize()
    pfes_truevsize = pfes.GetTrueVSize()

    un = mfem.Vector()
    un.SetSize(vfes_truevsize)
    un = 0.0
    un_next.SetSize(vfes_truevsize)
    un_next = 0.0
    unm1.SetSize(vfes_truevsize)
    unm1 = 0.0
    unm2.SetSize(vfes_truevsize)
    unm2 = 0.0
    fn.SetSize(vfes_truevsize)
     Nun.SetSize(vfes_truevsize)
     Nun = 0.0
     Nunm1.SetSize(vfes_truevsize)
     Nunm1 = 0.0
     Nunm2.SetSize(vfes_truevsize)
     Nunm2 = 0.0
     Fext.SetSize(vfes_truevsize)
     FText.SetSize(vfes_truevsize)
     Lext.SetSize(vfes_truevsize)
     resu.SetSize(vfes_truevsize)

    tmp1.SetSize(vfes_truevsize)

   pn.SetSize(pfes_truevsize);
   pn = 0.0;
   resp.SetSize(pfes_truevsize);
   resp = 0.0;
   FText_bdr.SetSize(pfes_truevsize);
   g_bdr.SetSize(pfes_truevsize);

   un_gf.SetSpace(vfes);
   un_gf = 0.0;
   un_next_gf.SetSpace(vfes);
   un_next_gf = 0.0;

   Lext_gf.SetSpace(vfes);
   curlu_gf.SetSpace(vfes);
   curlcurlu_gf.SetSpace(vfes);
   FText_gf.SetSpace(vfes);
   resu_gf.SetSpace(vfes);

   pn_gf.SetSpace(pfes);
   pn_gf = 0.0;
   resp_gf.SetSpace(pfes);

   cur_step = 0;

   PrintInfo();
    
    return mesh

