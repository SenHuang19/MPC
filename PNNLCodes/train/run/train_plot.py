import numpy as np
import pandas as pd
from itertools import repeat
import matplotlib.pyplot as plt

import sys


index=sys.argv[1]


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



i=0

for zone in zones:
    csp['z'+str(i)]=result[str(zone).upper()+':Zone Thermostat Cooling Setpoint Temperature [C](Each Call)']
    hsp['z'+str(i)]=result[str(zone).upper()+':Zone Thermostat Heating Setpoint Temperature [C](Each Call)']
    zt['z'+str(i)]=result[str(zone).upper()+':Zone Mean Air Temperature [C](TimeStep)']    
    m['z'+str(i)]=result['NODE '+str(termoutlet[i])+':System Node Mass Flow Rate [kg/s](TimeStep)']
    rh['z'+str(i)]=(result['NODE '+str(rhtinlet[i])+':System Node Temperature [C](TimeStep)']-result['NODE '+str(rhtoutlet[i])+':System Node Temperature [C](TimeStep)'])*result['NODE '+str(rhtoutlet[i])+':System Node Mass Flow Rate [kg/s](TimeStep)']*4200
    i=i+1    
      
       
#m=m.groupby(np.arange(len(m))//60).mean()

#rh=rh.groupby(np.arange(len(rh))//60).mean()


xlab=[]
xnum=[]
for j in range(5):
  for i in range(0,24,6):
    xlab.append(str(i)+':00')
    xnum.append((i*60)+24*j*60)

xlab.append('0:00')
xnum.append((5*24*60))




x=np.arange(len(m))


fig = plt.figure()

plt.subplot(2, 1, 1)
plt.plot(x,zt['z'+str(index)],label='Zone Temp',color='red')


plt.plot(x,csp['z'+str(index)],label='CSP',linestyle='-.',marker='x',markevery=120,color='gray')

plt.plot(x,hsp['z'+str(index)],label='HSP',linestyle='-.',marker='.', markevery=120,color='gray')
plt.xticks(xnum,[''],rotation=90)
plt.xlim(xnum[0],xnum[-1])
plt.ylim(18,28)
plt.legend(loc='best',ncol=3)
plt.ylabel('Temp [$^o$C]')

plt.subplot(2, 1, 2)


plt.plot(x,m['z'+str(index)])

plt.ylabel('Air Flow Rate [$kg/s$]')
# plt.xticks(xnum,[''],rotation=90)

# plt.subplot(3, 1, 3)


# plt.plot(x,rh['z0'],label='Zone 1')
# plt.plot(x,rh['z15'],label='Zone 2')
# plt.ylabel('Reheat [$W]')
plt.xticks(xnum,xlab,rotation=90)
plt.xlim(xnum[0],xnum[-1])
 

plt.savefig('train_'+str(index)+'.PNG',bbox_inches = 'tight',pad_inches = 0)
