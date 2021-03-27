'''

   vector_example.py

   demonstrate how to make mfem::Vector and assign values
   from numpy array

'''
import numpy as np
import mfem.ser as mfem
import traceback as tb

## Vector Creation
print("You can create mfem::Vector from numpy array")
vec = np.array([1,2,3,4,5.])
v = mfem.Vector(vec)
v.Print()

print("Please make sure that you are passing float")
vec = np.array([1,2,3,4,5]).astype(float)
v = mfem.Vector(vec)
v.Print()

print("This one does not work")
vec = np.array([1,2,3,4,5])
try:
   v = mfem.Vector(vec)
except:
   tb.print_exc()    

## Setting Vector using mfem::Vector or numpy array
vec = np.array([1,2,3,4,7.])
v.Assign(vec)
v.Print()

v2 = mfem.Vector(np.array([1,2,3,4,10.]))
v.Assign(v2)
v.Print()

print("This does not work (wrong type)")
try:
    vec = np.array([1,2,3,4,3])
    v.Assign(vec)
    v.Print()
except:
   tb.print_exc()    

print("Also, length must be the same")
try:
   vec = np.array([1,2,3,4.])
   v.Assign(vec)
   v.Print()      
except:
   tb.print_exc()


