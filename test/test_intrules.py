from __future__ import print_function
import os
import sys
import numpy as np

if len(sys.argv) > 1 and sys.argv[1] == '-p':   
    import mfem.par as mfem
    use_parallel = True
    from mfem.common.mpi_debug import nicePrint as print
    from mpi4py import MPI
    myid  = MPI.COMM_WORLD.rank
    
else:
    import mfem.ser as mfem
    use_parallel = False
    myid = 0
    
def run_test():
    a = mfem.IntegrationPointArray()
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.IntegrationPointArray(mt)
    assert a.Size() == 0, "array size is wrong"

    mt = mfem.MemoryType.HOST
    a = mfem.IntegrationPointArray(3, mt)
    assert a.Size() == 3, "array size is wrong"

    #If you do this. be carefull not to rules being GCed.
    #rules = mfem.IntegrationRules()
    #ir = rules.Get(2, 2)    
    ir = mfem.IntRules.Get(2, 2)

    answer = [(0.16666666666666666, 0.16666666666666666),
              (0.16666666666666666, 0.6666666666666667),
              (0.6666666666666667, 0.16666666666666666)]
    
    points = [(ir[i].x, ir[i].y) for i in range(ir.GetNPoints())]
    if points != answer:
        print(points)
        print(answer)
        assert False, "integration points coords are not collect" 
    
    pts = [ir[i] for i in range(ir.GetNPoints())]    
    a = mfem.IntegrationPointArray(pts)
    bbb = mfem.IntegrationPointArray((a.GetData(), 3))

    points2 = [(a[i].x, a[i].y) for i in range(ir.GetNPoints())]

    assert (points2 == points), "IntegrationPointArray coords does not agree (check 1)."

    irs = [mfem.IntRules.Get(i, 2) for i in range(mfem.Geometry.NumGeom)]
    
    rulearray = mfem.IntegrationRulePtrArray(irs)

    ir2 = rulearray[2]
    points3 = [(ir2[i].x, ir2[i].y) for i in range(ir2.GetNPoints())]
    assert (points3 == points), "IntegrationPointArray coords does not agree (check 2)."

    # check slice of IntegrationRulePtrArray
    rulearray2 = rulearray[:5]
    assert not rulearray2.OwnsData(), "subarray should not own data"
    ir2 = rulearray2[2]
    points3 = [(ir2[i].x, ir2[i].y) for i in range(ir2.GetNPoints())]
    assert (points3 == points), "IntegrationPointArray coords does not agree (check 3)."

    # create it from a pointer array
    data = rulearray2.GetData()    # this returns <Swig Object of type 'mfem::IntegrationRule **'>
    rulearray3 = mfem.IntegrationRulePtrArray((data, 5))
    ir2 = rulearray3[2]
    points3 = [(ir2[i].x, ir2[i].y) for i in range(ir2.GetNPoints())]
    assert (points3 == points), "IntegrationPointArray coords does not agree (check 3)."
    
    try:
        a = mfem.IntegrationPointArray([1,2,3])
    except:
        print("this one is supposed to fail")
        
if __name__=='__main__':
    run_test()
    
