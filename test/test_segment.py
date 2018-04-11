from __future__ import print_function
import os
import sys

def run_test():
    print("Test segment module")        
    seg = mfem.Segment()
    seg.SetVertices((1,2))
    seg.SetAttribute(1)
    assert  seg.GetVertices() == (1,2), "failed to set attribute (must be (1,2))"
    assert  seg.GetAttribute() == 1, "failed to set attribute (must be 1)"

    seg = mfem.Segment()
    seg.SetVertices((1,2))
    seg.SetAttribute(1)
    assert  seg.GetVertices() == (1,2), "failed to set attribute (must be (1,2))"
    assert  seg.GetAttribute() == 1, "failed to set attribute (must be 1)"

    seg = mfem.Segment([1, 2], 1)
    assert  seg.GetVertices() == (1,2), "failed to set attribute (must be (1,2))"
    assert  seg.GetAttribute() == 1, "failed to set attribute (must be 1)"

    seg = mfem.Segment((1, 2), 1)
    assert  seg.GetVertices() == (1,2), "failed to set attribute (must be (1,2))"
    assert  seg.GetAttribute() == 1, "failed to set attribute (must be 1)"

if __name__=='__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':   
        import mfem.par as mfem
    else:
        import mfem.ser as mfem
        
    run_test()

    



