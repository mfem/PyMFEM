from __future__ import print_function
import os
import sys

def run_test():
    print("Test point module")
    
    pt = mfem.Point(3)
    pt.SetAttribute(1)
    #print(pt.GetVertices(), pt.GetAttribute())
    assert  pt.GetVertices() == 3, "failed to set attribute (must be 3)"
    assert  pt.GetAttribute() == 1, "failed to set attribute (must be 1)"

if __name__=='__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':   
        import mfem.par as mfem
    else:
        import mfem.ser as mfem
    
    run_test()


