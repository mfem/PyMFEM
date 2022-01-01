'''

   run_tests.py

   run unit test cases

'''
import os
import sys
import subprocess as sp

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    
test_dir = os.path.dirname(os.path.abspath("__file__"))
sandbox = os.path.join(test_dir, 'sandbox')
os.makedirs(sandbox, exist_ok=True)

files = os.listdir(test_dir)

files = [x for x in files if x.endswith('.py')]
files = [x for x in files if not x.startswith('run_examples.py')]
files = [x for x in files if not x.startswith('run_tests.py')]
files = [x for x in files if not x.startswith('test_module.py')]
files = [x for x in files if not x.startswith('test_memory.py')]
                           

for t in  files:
   case = t[5:-3]

   d1 = os.path.join(sandbox, 's', case)
   d2 = os.path.join(sandbox, 'p', case)
   
   os.makedirs(d1, exist_ok=True)
   os.makedirs(d2, exist_ok=True)   

   os.chdir(d1)

   print("#### running (serial): " + t)

   comm = [sys.executable, "../../../"+t]
   
   p = sp.Popen(comm, stdout=sp.PIPE, stderr=sp.STDOUT, stdin=sp.PIPE)   
   lines, errs = p.communicate()
   print(lines.decode())   
   print(errs)

   os.chdir(d2)
   
   print("#### running (parallel): " + t)
   comm = [sys.executable,  "../../../"+t, "-p"]
   p = sp.Popen(comm, stdout=sp.PIPE, stderr=sp.STDOUT, stdin=sp.PIPE)   
   lines, errs = p.communicate()
   print(lines.decode())      
   print(errs)

