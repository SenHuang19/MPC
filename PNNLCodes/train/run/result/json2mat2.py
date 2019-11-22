import json
from scipy.io import savemat
import os

orig={}

for i in range(1,5):
           with open('fan'+str(i)+'.config') as handle:
                dictdump = json.loads(handle.read())
                temp={}
                for key in dictdump.keys():
#            print key.decode('utf8')
                       temp[key.encode("utf-8")]=dictdump[key]

           print temp
           orig['fan'+str(i)]=temp

savemat('fans.mat', orig, oned_as='row')
#print list(dictdump.keys()decode('utf8'))

with open('chiller.config') as handle:
                dictdump = json.loads(handle.read())
                temp={}
                for key in dictdump.keys():
#            print key.decode('utf8')
                       temp[key.encode("utf-8")]=dictdump[key]
                orig['chiller']=temp
savemat('chiller.mat', orig, oned_as='row')