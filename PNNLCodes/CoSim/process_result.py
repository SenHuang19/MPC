import numpy as np
import pandas as pd
from itertools import repeat
import matplotlib.pyplot as plt


price_e=np.loadtxt('paras/price_e.csv')

price_n=np.loadtxt('paras/price_g.csv')




ele=[float(x) for item in price_e for x in repeat(item, 3600/60)]

gas=[float(x) for item in price_n for x in repeat(item, 3600/60)]


result1=pd.read_csv('result/baseline/eplusout.csv')


result1['ele']=result1['CAV_BAS FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_BOT WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_MID WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['90.1-2010 WATERCOOLED  CENTRIFUGAL CHILLER 0 736TONS 0.6KW/TON:Chiller Electric Power [W](TimeStep)']

cost=0
fan=0
ch=0
boiler=0
for i in range(len(result1)):
              cost=cost+result1['ele'].iloc[i]*ele[i]
              fan=fan+result1['CAV_BAS FAN:Fan Electric Power [W](TimeStep)'].iloc[i]+result1['VAV_BOT WITH REHEAT FAN:Fan Electric Power [W](TimeStep)'].iloc[i]+result1['VAV_MID WITH REHEAT FAN:Fan Electric Power [W](TimeStep)'].iloc[i]+result1['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)'].iloc[i]
              ch=ch+result1['90.1-2010 WATERCOOLED  CENTRIFUGAL CHILLER 0 736TONS 0.6KW/TON:Chiller Electric Power [W](TimeStep)'].iloc[i]
              boiler=boiler+result1['HOT WATER LOOP BOILER 8806KBTU/HR 0.8 THERMAL EFF:Boiler Gas Rate [W](TimeStep)'].iloc[i]
              
baseline=cost/60./1000




opt=[]

actual=[]


sampling_times=[60,300,900,1800,3600]

for sampling_time in sampling_times:



     result1=pd.read_csv('result/3600_'+str(sampling_time)+'/eplusout.csv')


     result1['ele']=result1['CAV_BAS FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_BOT WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_MID WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)']+result1['90.1-2010 WATERCOOLED  CENTRIFUGAL CHILLER 0 736TONS 0.6KW/TON:Chiller Electric Power [W](TimeStep)']

     cost=0
     fan=0
     ch=0
     boiler=0
     for i in range(len(result1)):
              cost=cost+result1['ele'].iloc[i]*ele[i]
              fan=fan+result1['CAV_BAS FAN:Fan Electric Power [W](TimeStep)'].iloc[i]+result1['VAV_BOT WITH REHEAT FAN:Fan Electric Power [W](TimeStep)'].iloc[i]+result1['VAV_MID WITH REHEAT FAN:Fan Electric Power [W](TimeStep)'].iloc[i]+result1['VAV_TOP WITH REHEAT FAN:Fan Electric Power [W](TimeStep)'].iloc[i]
              ch=ch+result1['90.1-2010 WATERCOOLED  CENTRIFUGAL CHILLER 0 736TONS 0.6KW/TON:Chiller Electric Power [W](TimeStep)'].iloc[i]
              boiler=boiler+result1['HOT WATER LOOP BOILER 8806KBTU/HR 0.8 THERMAL EFF:Boiler Gas Rate [W](TimeStep)'].iloc[i]
              
     actual.append((baseline-cost/60./1000)/baseline*100)
     
     
     opt_cost=np.loadtxt('result/3600_'+str(sampling_time)+'/c.csv')
     
     opt.append((baseline-opt_cost.sum())/baseline*100)
     
x=np.arange(len(sampling_times))

plt.plot(sampling_times,opt,label='predicted',marker='o')

plt.plot(sampling_times,actual,label='actual',marker='x')

plt.legend()

plt.ylabel('Cost Saving [%]')

plt.xlabel('Model Discretization \n Interval [s]')

plt.savefig('result.PNG',bbox_inches = 'tight',pad_inches = 0)




