from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
'''

Test for PyMFEM

  This script runs C++ version and Python version examples
  and compare numbers in output.

  Usage
    usage: test.py [-h] [-np NP] [-verbose] [-parallel] [-serial]
               [-clean] [-ex EX] [-mfempdir '../../mfem']
               [-mfemsdir '../../mfem_ser']
    optional arguments:
      -h, --help          show this help message and exit
      -np NP              Number of processes for MPI
      -verbose            Verbose print out
      -parallel           Test parallel version
      -serial             Test serial version
      -clean              clean temporary files
      -ex EX              Test one example  python test.py -verbose
      -mfempdir           Parallel MFEM build directory (examples folder should exist)
      -mfemsdir           Serial MFEM build directory (examples folder should exist)

  Note:

  It compares numbers appears in the last 5 lines
  of output text, which seems sufficient.

  Number of precision is not the same between C++ and Python,
  the script only compare the first 3 significant digits.

  It ignores CPU seconds in output.

  

'''
import os
import sys
import subprocess as sp
import time
import shutil
import re

skip_test = [''] 
#skip_test = ['ex15p']  #skip this since it takes too long on my Mac..;D

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

## special handling of text (needs to be set manually)
ignore_txt_def = ['seconds', 'residual']    
ignore_txt_dict = {'ex15p':['seconds', 'residual',],
                   'ex10':['iteration', 'seconds', 'residual',],
                   'ex10p':['iteration', 'seconds', 'residual',],
                   'ex13p':['iteration', 'seconds', 'residual',]}


def run_file(command, num = 5):
    t1 = time.time()
    p = sp.Popen(command,stdout=sp.PIPE, stderr=sp.STDOUT, stdin=sp.PIPE)
    p.stdin.write(b'q\n')
    lines, errs = p.communicate()
    lines = lines.decode('utf-8').split('\n')
    t2 = time.time()
    sys.stdout.flush()
    lines = [l for l in lines if len(l.strip()) != 0]
    return t2-t1, lines[-num:]

def pad_lines(lines):
    length = max([len(l.strip()) for l in lines])
    return  [l.strip() + " "*(length-len(l.strip())) for l in lines]

def extract_numbers(line, case):
    nums = []
    line = line.strip()

    igt = ignore_txt_dict.get(case, ignore_txt_def)
    for x in igt: 
        if line.lower().find(x) != -1: return []
    for x in line.split():
        try:
            if x.endswith(','): x = x[:-1]
            if x.endswith('.'): x = x[:-1]            
            x = float(x)
            nums.append("{:.3e}".format(x))
        except:
            pass
    return nums

def compare_results(lines1, lines2, case, verbose=True):
    lines1 = pad_lines(lines1)
    lines2 = pad_lines(lines2)

    compare = []
    for l1, l2 in zip(lines1, lines2):
        nums1 = extract_numbers(l1, case)
        nums2 = extract_numbers(l2, case)
        flag = nums1 == nums2
        if len(nums1) == 0 and len(nums2) == 0: flag = 'ignored'
        if verbose:print(l1 + "\t" + l2 + " : " + str(flag))
        if flag == 'ignored': continue
        compare.append(flag)
    return all(compare)
    
def run_test(mfem_exes, pymfem_exes, serial=True, np=2, verbose=False):
    print("mfem examples from : " + os.path.dirname(mfem_exes[0]))
    print("PyMFEM examples from : " + os.path.dirname(pymfem_exes[0]))
    print(',  '.join([os.path.basename(x) for x in mfem_exes]))    

    if not serial:
        comm_mpi = ["mpirun", "-np", str(np)]
    else:
        comm_mpi = []

    results = []; fails = []

    for e1, e2 in zip(mfem_exes, pymfem_exes):
        case =  os.path.basename(e1)
        print("Running : " + case)
        comm = comm_mpi + [e1]
        t1, l1 = run_file(comm, num = 5)

        print("Running : " + os.path.basename(e2))

        # note -u is unbuffered option
        comm = comm_mpi + [sys.executable,  "-u", e2]
        t2, l2 = run_file(comm, num = 5)

        flag = compare_results(l1, l2, case, verbose=verbose)
        if flag:
            print(os.path.basename(e1) + "  " +  bcolors.OKGREEN + "PASS" + bcolors.ENDC)
        else:
            print(os.path.basename(e1) + "  " +  bcolors.FAIL + "FAIL" + bcolors.ENDC)
            fails.append(os.path.basename(e1))
        results.append(flag)
        print("Time spent for execution : (C) " + str(t1) + " (Python) " + str(t2) + " " + str((t2/t1-1)*100) + "%")
        
    return results, fails
        
def find_mfem_examples(dir, serial = True, example='all'):
    if dir == '': return []
    example_dir = os.path.join(dir, 'examples')
    names = [x for x in os.listdir(example_dir)
             if x.startswith('ex') and not x.endswith('.cpp')]

    if serial:
        names = [x for x in names if not x.endswith('p')]
    else:    
        names = [x for x in names if x.endswith('p')]

    if example != 'all':
        names = [n for n in names if n == example]
    names = [n for n in names if not n in skip_test]    

    nums = [int(''.join(re.findall(r'\d+', n))) for n in names]
    names = [os.path.join(example_dir, x) for x in names]    
    names = [x for x in names if os.access(x, os.X_OK)]

    pairs = [(num, name) for num, name in zip(nums, names)]
    names = [name  for num, name in sorted(pairs)]

    return names

def print_help(self):
    print("python test.py -np -v")

if __name__=="__main__":
    script_path = sys.path[0]  # the location of this file
    
    from mfem.common.arg_parser import ArgParser
    parser = ArgParser(description='test')
    parser.add_argument('-np', default = 2,
                        action = 'store', type = int,
                    help='Number of processes for MPI')
    parser.add_argument('-verbose', 
                        action = 'store_true', default = False,
                        help='Verbose print out')
    parser.add_argument('-parallel', 
                        action = 'store_true', 
                        help='Test parallel version')
    parser.add_argument('-serial', 
                        action = 'store_true', 
                        help='Test serial version')
    parser.add_argument('-clean', 
                        action = 'store_true', default = False,
                        help='Test serial version')
    parser.add_argument('-ex', 
                        action = 'store', default = "all",
                        help='Test one example')
    parser.add_argument('-mfempdir', 
                        action = 'store', default = "../../mfem",
                        help='mfem (parallel) directory')
    parser.add_argument('-mfemsdir', 
                        action = 'store', default = "../../mfem_ser",
                        help='mfem (serial) directory')
    
    
    
    args = parser.parse_args()
    parser.print_options(args)

    np = args.np
    verbose = args.verbose
    testp = bool(args.parallel)
    tests = bool(args.serial)
    clean = args.clean
    example = args.ex
    mfempdir = args.mfempdir
    mfemsdir = args.mfemsdir    


    if clean:
        od = os.getcwd()
        files = os.listdir(od)
        for x in files:
           if x.endswith(".py"): continue
           if os.path.isdir(x):
               shutil.rmtree(x)
           else:
               os.remove(x)
        sys.exit()
    import mfem

    dir1 = os.path.dirname(script_path)
    dir2 = os.path.dirname(os.path.dirname(mfem.__file__))

    print(bcolors.FAIL+'PyMFEM module loaded from :' + dir2 + bcolors.ENDC)
    print(bcolors.FAIL+'Searching Python Example in :' + dir1 + bcolors.ENDC)

    exdir = os.path.join(dir1, "examples")
    
    mfem_dir = os.path.abspath(mfempdir)
    mfemser_dir = os.path.abspath(mfemsdir)

    if testp:
       mfemp_exes = find_mfem_examples(mfem_dir, serial = False, example = example)
       pymfemp_exes = [os.path.join(exdir, os.path.basename(x))+".py"  for x in mfemp_exes]
    if tests:
       mfems_exes = find_mfem_examples(mfemser_dir, serial = True,  example = example)    
       pymfems_exes = [os.path.join(exdir, os.path.basename(x))+".py"  for x in mfems_exes]


    if tests and mfems_exes:
        print(bcolors.BOLD+"Running Serial Version"+bcolors.ENDC)
        results, fails = run_test(mfems_exes, pymfems_exes, serial = True, verbose=verbose)
    else:
        results = []; fails=[]        
    if testp and len(mfemp_exes) != 0:
        print(bcolors.BOLD+"Running Parallel Version"+bcolors.ENDC)
        resultp, failp = run_test(mfemp_exes, pymfemp_exes, serial = False, np = np, verbose=verbose)
    else:
        resultp = []; failp=[]


    if tests and len(results) != 0:
        if all(results):
            print("Serial Test \t" + bcolors.OKGREEN + "ALL PASS" + bcolors.ENDC)
        else:
            print("Serial Test \t" + bcolors.FAIL + " ".join(fails) + bcolors.ENDC)

    if testp and len(resultp) != 0:
        if all(resultp):
            print("Parallel Test \t" + bcolors.OKGREEN + "ALL PASS" + bcolors.ENDC)
        else:
            print("Parallel Test \t" + bcolors.FAIL + " ".join(failp) + bcolors.ENDC)        
        
    
             
    
    
    
    
