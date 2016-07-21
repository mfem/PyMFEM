#
#  Helper script to create setup_local.py
#
import os
import sys

if sys.version_info < (3, 0):
    def execute(text, l = {}, g = {}):
        exec text in l, g
else:
    def execute(text, l = {}, g = {}):
        exec(text, l, g)

fid = open('setup_local.py', 'w')
fid.write("#  setup_local.py \n")
fid.write("#  generated from write_setup_local.py\n")
fid.write("#  do not edit this directly\n")
fid.write("#  instead use Make setup_local.py\n")
for key in os.environ:
    try:
        text = key.lower() + ' = "' + os.environ[key] + '"'
        execute(text, {}, {})
        fid.write(text+"\n")
    except:
        pass
fid.close()    
