%define ENUM_CLASS_WRAP(P1, P2)
%pythoncode %{
 
import sys
__thismodule = sys.modules[__name__]

from enum import IntEnum
  
__enum = P1             # your enumeration
__scope_name = "P2" + "_"   # scope's name 
__scope_length = len(__scope_name) # length
__values = {}  
  
for name in dir(__enum):
    if name.find(__scope_name) == 0:
        __values[name[__scope_length:]] = getattr(__enum, name)
        delattr( __thismodule, __scope_name + name[__scope_length:])

P2 = IntEnum("P2", __values)
  
del name, __enum, __scope_name, __scope_length, __thismodule, __values
%}
%enddef
