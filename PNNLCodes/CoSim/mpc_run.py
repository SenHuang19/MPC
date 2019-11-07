import socket
import matlab.engine
import numpy as np
from itertools import repeat
import sys
import pandas as pd
import random as rd

class socket_server:

    def __init__(self):     
          self.sock=socket.socket()
          self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
          host=socket.gethostname()
          port=47569
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

write_port_file(47569,'127.0.0.1')

vers = 2
flag = 0

### read price info

price_e=np.loadtxt('price_e.csv')

price_n=np.loadtxt('price_g.csv')

control_horizon_length=300

samp_time=300


ele=[float(x) for item in price_e for x in repeat(item, 3600/samp_time)]

gas=[float(x) for item in price_n for x in repeat(item, 3600/samp_time)]


print len(ele)

### read prediction

prediction=pd.read_csv('input.csv' )


prediction=prediction.groupby(np.arange(len(prediction))//int(samp_time/60)).mean()

tout=prediction['tout'].values.tolist()

print len(tout)


          
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
              if index%(int(control_horizon_length/60))==0:

                       start_t=index/int(control_horizon_length/60)  

                       end_t=start_t+1440/int(control_horizon_length/60)    
                       
                       print index
             
                       print start_t
                       
                       print end_t
                                              
                       
                       ele_m=matlab.double(ele[start_t:end_t])

                       gas_m=matlab.double(gas[start_t:end_t])
         
                       tini=[]
                       for i in range(7,7+16):
                                  tini.append(float(arry[i]))                                       
                       tini_m=matlab.double([tini])         
         
         
                       qin=[]

                       for i in range(1,6):
                            qin.append(prediction['top'+str(i)].values.tolist()[start_t:end_t])

                       for i in range(1,6):
                            qin.append(prediction['mid'+str(i)].values.tolist()[start_t:end_t])

                       for i in range(1,6):
                            qin.append(prediction['bot'+str(i)].values.tolist()[start_t:end_t])
    
                       qin.append(prediction['basement'].values.tolist()[start_t:end_t])
         
                       qin_m=matlab.double(qin)

                       tout_m=matlab.double(tout[start_t:end_t])             
              
#                          print index
                       set=[]
                       tset= eng.func_EDC_CoSim_test(samp_time/60.,ele_m,gas_m,tout_m,qin_m,tini_m)
#                       print tset

                      
              for i in range(16):
                   mssg = mssg + ' ' + str(tset[i][0])+ ' ' + str(tset[i][0])
              mssg = mssg+ ' ' + str(0.3) + ' ' + str(0.3)+ ' ' +str(0.3)+ ' ' +str(0.3)+ ' ' +str(2500)+ ' ' +str(0)

              mssg =  mssg+'\n'
              index=index+1
              conn.send(mssg)

	 

	 

