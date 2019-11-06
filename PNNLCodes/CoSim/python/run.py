import os
import subprocess
import time
from shutil import copyfile
import sys
import json,collections


modelPath=''

os.environ['BCVTB_HOME']='bcvtb'


weatherPath='USA_MD_Baltimore-Washington.Intl.AP.724060_TMY3.epw'

modelDir='B_Office_20181228.idf'


cmdStr = "C:\EnergyPlusV8-9-0\energyplus -w \"%s\" -r \"%s\"" % (weatherPath, modelDir)

simulation = subprocess.Popen(cmdStr, shell=True)
sock=subprocess.Popen('python '+'mpc_run.py', shell=True)


sock.wait()
#sock.terminate()
