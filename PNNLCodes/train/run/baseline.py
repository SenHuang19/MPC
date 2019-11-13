import socket
import matlab.engine
import numpy as np
from itertools import repeat
import sys
import pandas as pd
import random as rd
import math

class socket_server:

    def __init__(self):     
          self.sock=socket.socket()
          self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
          host=socket.gethostname()
          port=48000
          self.sock.bind(("127.0.0.1",port))
          self.sock.listen(10)
		  
def data_parse(data): 
    data=data.replace('[','')     
    data=data.replace(']','')  	
    data=data.split(',')

    for i in range(len(data)):
	              data[i]=float(data[i])
         
    return data




def write_port_file(port,host):
        fh = open('socket.cfg', "w+")
        fh.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        fh.write('<BCVTB-client>\n')
        fh.write('  <ipc>\n')
        fh.write('    <socket port="%r" hostname="%s"/>\n' % (port, host))
        fh.write('  </ipc>\n')
        fh.write('</BCVTB-client>')
        fh.close()
		  	 		  

### start the matlab engine (optional)

eng = matlab.engine.start_matlab()

### set up the interface

server=socket_server()

write_port_file(48000,'127.0.0.1')

vers = 2
flag = 0


          
server.sock.listen(10)


conn,addr=server.sock.accept()
index=0


#print len(samp_time)




while 1:



         
#         print len(qin[1][index*int(control_horizon_length/60):index*int(control_horizon_length/60)+1440/int(control_horizon_length/60)])

         

### data received from energyplus
         
         data = conn.recv(10240)
		 
         data = data.rstrip()

         arry = data.split()
         flagt = float(arry[1])
         if flagt==1:
                 conn.close()
                 sys.exit()
         if len(arry)>6:
              time=float(arry[5])
              mssg = '%r %r %r 0 0 %r' % (vers, flag, 38, time) 

                
#              print mssg

              if index<=1439:
                   tset=22
              else:
                   tset=22              
              for i in range(16):
                   mssg = mssg + ' ' + str(20)+ ' ' + str(22+2*math.sin(index/1440.*2*3.14))
              mssg = mssg+ ' ' + str(0.1) + ' ' + str(0.1)+ ' ' +str(0.1)+ ' ' +str(0.1)+ ' ' +str(0)+ ' ' +str(0)

              mssg =  mssg+'\n'
              index=index+1
              conn.send(mssg)

	 

	 

