import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
import sys

import json
import collections
from collections import OrderedDict
import itertools as it
import math
import pandas as pd

from scipy.optimize import minimize 
from scipy.optimize import nnls 
from sklearn.linear_model import Lasso

import json
import collections
from collections import OrderedDict

regr = linear_model.LinearRegression(fit_intercept = False)

lin = Lasso(alpha=0.0001,precompute=True,max_iter=1000,
            positive=True, random_state=9999, selection='random')

s=sys.argv[1]
filename='zone'+sys.argv[2]

raw=pd.read_csv('raw/'+s+'_'+filename+'.csv')
  
#train the model
tab=raw[:1440*4]

#control horizon length
num=5

print len(tab)

tab=tab.groupby(np.arange(len(tab))//num).mean()

print len(tab)

# the equation we have in mind is 
# c(T^k+1-T^k)=(Ta^k-T^k)/R+mdtC^air+RH+I+W
# T^k+1=Ta^k/c+(1-1/RC)T^k+mdtC^air/c+RH/c+I/c+W/c
# T^k+1=a_1Ta^k+a_2T^k+a_3m+a_4RH+a_5I+a_0W

x=zip(tab['tout'][1:],tab['t1'][:-1],tab['m'][1:],tab['rh'][1:],tab['i'][1:])
y=tab['t1'][1:]

regr.fit(x, y)

x=zip(tab['tout'][1:],tab['t1'][:-1],-tab['m'][1:],tab['rh'][1:],tab['i'][1:])

lin.fit(x, y)

print lin.coef_

a1=lin.coef_[0]

a2=lin.coef_[1]

a3=-lin.coef_[2]

a4=lin.coef_[3]

a5=lin.coef_[4]

a0=lin.intercept_ 

y_pred=[]

tab=raw[1440*4:1440*5]

tab.to_csv('temp.csv')

tab=pd.read_csv('temp.csv')

tab=tab.groupby(np.arange(len(tab))//num).mean()
print len(tab)

t_pre=tab['t1'].iloc[0]

for i in range(1,len(tab)):
#    t_pre=tab['t1'].iloc[i-1]
    t_pre=a1*tab['tout'][i]+a2*t_pre+a3*tab['m'].iloc[i]+a4*tab['rh'].iloc[i]+a5*tab['i'].iloc[i]+a0
    y_pred.append(t_pre)

y=[]
for i in range(1,len(tab)):
      y.append(tab['t1'].iloc[i])



tab3=pd.DataFrame()
tab3['y_pred']=y_pred

tab3['real']=y

xx=np.arange(len(tab3))
plt.plot(xx,tab3['y_pred'], label='prediction',marker='x',markevery=20)


plt.plot(xx,tab3['real'], label='real')
plt.legend()
#plt.xticks(xtic,xlab,rotation=45)
plt.xlim(xx[0],xx[-1])
plt.savefig('result/'+s+'_'+filename+'testlidation_result.png',bbox_inches = 'tight',pad_inches = 0.1)
plt.show()

coefs=collections.OrderedDict()
coefs['a0']=a0
coefs['a1']=a1
coefs['a2']=a2
coefs['a3']=a3
coefs['a4']=a4
coefs['a5']=a5

with open('result/'+filename+'_'+str(s)+'_result','w') as outfile:
      json.dump(coefs,outfile, sort_keys=True, indent=2) 


