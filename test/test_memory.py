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

import resource

def format_memory_usage(point="memory usage"):
    usage=resource.getrusage(resource.RUSAGE_SELF)
    txt = "%s: usertime=%s systime=%s mem="+ bcolors.FAIL + "%s" + bcolors.ENDC + " mb"
    return txt%(point,usage[0],usage[1],
               (usage[2]*resource.getpagesize())/1000000.0 )


if __name__ == '__main__':
   import ex1p

   for i in range(30):
       ex1p.run()
       print(format_memory_usage())
