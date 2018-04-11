from __future__ import print_function

import os
import sys
import time
import subprocess as sp

modules = ["point", "segment", "mesh"]


def run_file(command):
    print('command',command)
    t1 = time.time()    
    p = sp.Popen(command,stdout=sp.PIPE, stderr=sp.STDOUT, stdin=sp.PIPE)
    lines = p.stdout.readlines()
    t2 = time.time()
    sys.stdout.flush()
    print(lines)
    return t2-t1, lines

def run_test(serial = True, np=2):
    file = ["test_"+x+".py" for x in modules]
    
    if not serial:
        comm_mpi = ["mpirun", "-np", str(np)]
    else:
        comm_mpi = []

    for f in file:
        comm = comm_mpi + [sys.executable,  f]
        if not serial:
            comm = comm.append("-p")


        run_file(comm)

if __name__=='__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '-p':   
        run_test(serial = False)
    else:
        run_test(serial = True)

