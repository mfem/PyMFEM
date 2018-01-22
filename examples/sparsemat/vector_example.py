import numpy as np
import mfem.ser as mfem


vec = np.array([1,2,3,4,5])
try:
   v = mfem.Vector(vec)
except:
   print("can not give int")
   
print("pass float")
vec = np.array([1,2,3,4,5]).astype(float)
v = mfem.Vector(vec)
v.Print()


print("or this way")
vec = np.array([1,2,3,4,5]).astype(float)
v = mfem.Vector(vec)
v.Print()

vec = np.array([1,2,3,4,7])
v.Assgin(vec)

