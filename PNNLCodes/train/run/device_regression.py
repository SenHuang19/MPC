import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
import sys


import itertools as it
import math
import pandas as pd

from scipy.optimize import minimize 
from scipy.optimize import nnls 
from sklearn.linear_model import Lasso
from scipy.optimize import lsq_linear
import json
import collections
from collections import OrderedDict

from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error


regr = linear_model.LinearRegression(fit_intercept = True)

lin = Lasso(alpha=0.0001,precompute=True,max_iter=1000,
            positive=True, random_state=9999, selection='random')
            
oa= pd.read_csv('result/oa.csv')

m = pd.read_csv('result/m.csv')

rh = pd.read_csv('result/rh.csv')

zt = pd.read_csv('result/zt.csv')


rt=[22*0.9]*len(oa)

st=[(55-32)/9.*5]*len(oa)
t1=oa['oa']*0.1


dt=t1+rt-st



device = pd.read_csv('result/device.csv')


m_tot1=m['z'+str(0)]
for i in range(1,5):
      m_tot1=m_tot1+m['z'+str(i)]

x=zip(m_tot1,m_tot1*m_tot1)

y=device['fan1']/1000
regr.fit(x, y)

y2=regr.predict(x)
plt.clf()


plt.scatter(y,y2,marker='x')

plt.savefig('result/fan1.png',bbox_inches = 'tight',pad_inches = 0.1)





coefs=collections.OrderedDict()
coefs['c0']=regr.intercept_ 
coefs['c1']=regr.coef_[0]
coefs['c2']=regr.coef_[1]

with open('result/fan1.config','w') as outfile:
          json.dump(coefs,outfile, sort_keys=True, indent=2) 


m_tot2=m['z'+str(5)]
for i in range(6,10):
      m_tot2=m_tot2+m['z'+str(i)]

x=zip(m_tot2,m_tot2*m_tot2)

y=device['fan2']/1000

regr.fit(x, y)

y2=regr.predict(x)


plt.clf()


plt.scatter(y,y2,marker='x')

plt.savefig('result/fan2.png',bbox_inches = 'tight',pad_inches = 0.1)

coefs=collections.OrderedDict()
coefs['c0']=regr.intercept_ 
coefs['c1']=regr.coef_[0]
coefs['c2']=regr.coef_[1]

with open('result/fan2.config','w') as outfile:
          json.dump(coefs,outfile, sort_keys=True, indent=2) 

m_tot3=m['z'+str(10)]
for i in range(11,15):
      m_tot3=m_tot3+m['z'+str(i)]

x=zip(m_tot3,m_tot3*m_tot3)

y=device['fan3']/1000

regr.fit(x, y)
y2=regr.predict(x)

plt.clf()



plt.savefig('result/fan3.png',bbox_inches = 'tight',pad_inches = 0.1)

coefs=collections.OrderedDict()
coefs['c0']=regr.intercept_ 
coefs['c1']=regr.coef_[0]
coefs['c2']=regr.coef_[1]

with open('result/fan3.config','w') as outfile:
          json.dump(coefs,outfile, sort_keys=True, indent=2) 



m_tot4=m['z'+str(15)]

x=zip(m_tot4)

y=device['fan4']/1000

regr.fit(x, y)
y2=regr.predict(x)


plt.clf()



plt.scatter(y,y2,marker='x')

plt.savefig('result/fan4.png',bbox_inches = 'tight',pad_inches = 0.1)

coefs=collections.OrderedDict()
coefs['c0']=regr.intercept_ 
coefs['c1']=regr.coef_[0]
coefs['c2']=0
with open('result/fan4.config','w') as outfile:
          json.dump(coefs,outfile, sort_keys=True, indent=2) 



qch=(m_tot1+m_tot2+m_tot3+m_tot4)*dt*1000.



x=zip(oa['oa'],qch)

y=device['chiller']/1000.

regr.fit(x, y)

print regr.coef_

print regr.intercept_


y2=regr.predict(x)


plt.clf()

plt.scatter(y,y2,marker='x')

plt.savefig('result/chiller.png',bbox_inches = 'tight',pad_inches = 0.1)


coefs=collections.OrderedDict()
coefs['d0']=regr.intercept_ 
coefs['d1']=regr.coef_[0]
coefs['d2']=regr.coef_[1]

with open('result/chiller.config','w') as outfile:
          json.dump(coefs,outfile, sort_keys=True, indent=2) 

rh_tot=rh['z'+str(10)]
for i in range(0,15):
      rh_tot=rh_tot+rh['z'+str(i)]

print rh_tot














