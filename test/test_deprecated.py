import sys
import warnings
warnings.simplefilter('default')

def test(parallel=False):
    if parallel:
        import mfem.par as mfem
    else:
        import mfem.ser as mfem


    pt = mfem.intp()
    
    tri = mfem.Triangle()
    tri.GetNFaces(pt)
    
if __name__ == "__main__":
    if "-p" in sys.argv:
        test(True)
    else:
        test(False)
        
