from __future__ import print_function
import os
import sys


def clean(dryrun = True):
    cdir = os.path.dirname(os.path.abspath(__file__))
    pdir = os.path.join(cdir, 'mfem', 'par')

    py_file = [x for x in os.listdir(pdir) if x.endswith('.py')]
    py_file = [x for x in py_file if not x.startswith('.')]


    par_file = [x for x in py_file if x.startswith('p')]

    ser_file = [x for x in py_file if not x.startswith('p')]

    names = ['point']
    for n in names:
        par_file.remove('point.py')
        ser_file.append('point.py')    

    par_import = ['import '+x[:-3] for x in par_file]

    for f in ser_file:
        fid = open(os.path.join(pdir, f), 'r')
        lines = fid.readlines()    
        print(f, [l.strip() for l in lines  if l.strip() in par_import])
        fid.close()
    if dryrun:
        return
    for f in ser_file:
        fid = open(os.path.join(pdir, f), 'r')
        lines = fid.readlines()
        fid.close()    
        lines = ['' if l.strip() in par_import else l  for l in lines]
        fid = open(os.path.join(pdir, f), 'w')
        fid.write(''.join(lines))
        fid.close()

if __name__=='__main__':
    dryrun = True
    if len(sys.argv) > 1:
        if sys.argv[1] == '-x': dryrun = False
    else:
        print('this is dry run. Use clean_import.py -x to perform cleaning')
    clean(dryrun)
    
