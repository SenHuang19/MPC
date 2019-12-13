import numpy as np
import pandas as pd
from itertools import repeat
import matplotlib.pyplot as plt
import sys

num=int(sys.argv[1])


result=pd.read_csv('eplusout.csv')


zones=['Perimeter_top_ZN_1 ZN','Perimeter_top_ZN_2 ZN','Perimeter_top_ZN_3 ZN','Perimeter_top_ZN_4 ZN','Core_top ZN','Perimeter_mid_ZN_1 ZN','Perimeter_mid_ZN_2 ZN','Perimeter_mid_ZN_3 ZN','Perimeter_mid_ZN_4 ZN','Core_mid ZN','Perimeter_bot_ZN_1 ZN','Perimeter_bot_ZN_2 ZN','Perimeter_bot_ZN_3 ZN','Perimeter_bot_ZN_4 ZN','Core_bottom ZN','Basement ZN']

termoutlet=[151,168,173,178,183,112,129,134,139,144,73,90,95,100,105,190]

rhtoutlet=[162,167,172,177,182,143,128,133,138,143,84,89,94,99,104,193]

rhtinlet=[161,165,171,176,181,142,126,132,137,142,83,87,93,98,103,192]




csp=pd.DataFrame()

hsp=pd.DataFrame()

m=pd.DataFrame()

zt=pd.DataFrame()

rh=pd.DataFrame()

q=pd.DataFrame()

oa=pd.DataFrame()
oa['oa']=result['Environment:Site Outdoor Air Drybulb Temperature [C](TimeStep)']

device=pd.DataFrame()

device['fan4']=result['CAV_BAS FAN:Fan Electric Power [W](TimeStep)']

device['fan3']=result['VAV_BOT WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']

device['fan2']=result['VAV_MID WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']

device['fan1']=result['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']

device['chiller']=result['90.1-2010 WATERCOOLED  CENTRIFUGAL CHILLER 0 736TONS 0.6KW/TON:Chiller Electric Power [W](TimeStep)']


i=0
min=[]
for zone in zones:
    if i>4 and i<10:
	    mut=10
    else:
	    mut=1	
    csp['z'+str(i)]=result[str(zone).upper()+':Zone Thermostat Cooling Setpoint Temperature [C](Each Call)']
            
    q['z'+str(i)]=result[str(zone).upper()+':Zone Air Heat Balance Internal Convective Heat Gain Rate [W](TimeStep)']
    hsp['z'+str(i)]=result[str(zone).upper()+':Zone Thermostat Heating Setpoint Temperature [C](Each Call)']
    zt['z'+str(i)]=result[str(zone).upper()+':Zone Mean Air Temperature [C](TimeStep)']    
    m['z'+str(i)]=result['NODE '+str(termoutlet[i])+':System Node Mass Flow Rate [kg/s](TimeStep)']/mut
    min.append(m['z'+str(i)].min())
    rh['z'+str(i)]=(result['NODE '+str(rhtinlet[i])+':System Node Temperature [C](TimeStep)']-result['NODE '+str(rhtoutlet[i])+':System Node Temperature [C](TimeStep)'])*result['NODE '+str(rhtoutlet[i])+':System Node Mass Flow Rate [kg/s](TimeStep)']*4200/mut
    i=i+1    
      
print min 
f=open('min.csv','w')
for i in range(len(min)):
    f.writelines(str(min[i])+'\n')
f.close()

       
m=m.groupby(np.arange(len(m))//num).mean()

rh=rh.groupby(np.arange(len(rh))//num).mean()

q=q.groupby(np.arange(len(q))//num).mean()

zt=zt.groupby(np.arange(len(zt))//num).mean()

oa=oa.groupby(np.arange(len(oa))//num).mean()

device=device.groupby(np.arange(len(device))//num).mean()

m.to_csv('result/m.csv')

rh.to_csv('result/rh.csv')

q.to_csv('result/q.csv')

zt.to_csv('result/zt.csv')

oa.to_csv('result/oa.csv')

device.to_csv('result/device.csv')
