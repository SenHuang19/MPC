import numpy as np
import pandas as pd
from itertools import repeat


price_e=np.loadtxt('paras/price_e.csv')

price_n=np.loadtxt('paras/price_g.csv')

ele=[float(x) for item in price_e for x in repeat(item, 3600/60)]

gas=[float(x) for item in price_n for x in repeat(item, 3600/60)]

print str('test')+' '+str('fan')+' '+str('ch')+' '+str('boiler')+' '+str('cost')


result1=pd.read_csv('eplusout.csv')


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
print str('opt')+' '+str(fan)+' '+str(ch)+' '+str(boiler)+' '+str(cost/60./1000)


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
print str('base')+' '+str(fan)+' '+str(ch)+' '+str(boiler)+' '+str(cost/60./1000)



