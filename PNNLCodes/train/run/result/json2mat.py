import json
from scipy.io import savemat
import os

orig={}

for i in range(16):
           with open('z'+str(i)+'.config') as handle:
                dictdump = json.loads(handle.read())
                temp={}
                for key in dictdump.keys():
#            print key.decode('utf8')
                       temp[key.encode("utf-8")]=dictdump[key]

           print temp
           orig['z'+str(i)]=temp
print orig
savemat('zones.mat', orig, oned_as='row')
#print list(dictdump.keys()decode('utf8'))