import socket
import time
import matlab.engine
import numpy as np
from itertools import repeat
import random as rd
import sys

class socket_server:

    def __init__(self):     
          self.sock=socket.socket()
          self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
          host=socket.gethostname()
          port=5000
          self.sock.bind(("127.0.0.1",5500))
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
		  
def callfun(tz,to,ts,reg):

   Tz=[]
   for i in range(len(tz)):
         Tz.append([tz[i]])
	  
	  
   Tz=matlab.double(Tz)
   


   To=[to,to,to,to,to]

   To=matlab.double([To])

   Ts=[]
   for i in range(len(tz)):
         Ts.append([ts[i],ts[i],ts[i],ts[i],ts[i]])			   

 
   Ts=matlab.double(Ts)	  

   reg=matlab.double([reg])
   
   e=[[1],[2],[3],[4],[5],[6],[7],[8],[9]]
   
   

   e=eng.mdot_calc_z05_10min(Tz,To,Ts,reg)
   

     
   return e		  
		  
		  


### start the matlab engine (optional)

#eng = matlab.engine.start_matlab()

server=socket_server()

### user define if the model is baseline 0 or 1

ele_raw=np.loadtxt('prices.csv')

control_horizon_length=300

ele=[x for item in ele_raw for x in repeat(item, 3600/control_horizon_length)]

#print reg
	
server=socket_server()
write_port_file(5500,'127.0.0.1')

vers = 2
flag = 0


          
server.sock.listen(10)

conn,addr=server.sock.accept()
index=0
while 1:


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

         
              print mssg
              if index%(int(control_horizon_length/60))==0:
                          print index
                          set=[]
                          for i in range(16):
                                     set.append(round(24+np.sin(index/288.0*2*3.14*(rd.uniform(-0.5,0.5)+1))*1,2))                          
              for i in range(16):
                   mssg = mssg + ' ' + str(set[i])+ ' ' + str(set[i])
              mssg = mssg+ ' ' + str(0.3) + ' ' + str(0.3)+ ' ' +str(0.3)+ ' ' +str(0.3)+ ' ' +str(2500)+ ' ' +str(0)

              mssg =  mssg+'\n'
              index=index+1
              conn.send(mssg)

	 

	 

