import json
from scipy.io import savemat
import os

orig={}

file_list=os.listdir(r"result")
#print (file_list)

for file in file_list:
    if file.find('result')!=-1 and file.find('png')==-1 and file.find('csv')==-1:
           print file
           with open('result/'+file) as handle:
                dictdump = json.loads(handle.read())
                temp={}
                for key in dictdump.keys():
#            print key.decode('utf8')
                       temp[key.encode("utf-8")]=dictdump[key]
#            print dictdump[key.decode('utf8')]

           print temp
           orig[file.replace('-','_')]=temp
print orig
savemat('mid_floor.mat', orig, oned_as='row')
#print list(dictdump.keys()decode('utf8'))