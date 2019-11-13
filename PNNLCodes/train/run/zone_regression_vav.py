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

regr = linear_model.LinearRegression(fit_intercept = True)

lin = Lasso(alpha=0.0001,precompute=True,max_iter=1000,
            positive=True, random_state=9999, selection='random')


z_no='z'+sys.argv[1]

m = pd.read_csv('result/m.csv')

rh = pd.read_csv('result/rh.csv')

q = pd.read_csv('result/q.csv')

zt = pd.read_csv('result/zt.csv')

oa= pd.read_csv('result/oa.csv')

# the equation we have in mind is 
# c(T^k+1-T^k)=(Ta^k-T^k)/R+mdtC^air+RH+I+W
# T^k+1=Ta^k/c+(1-1/RC)T^k+mdtC^air/c+RH/c+I/c+W/c
# T^k+1=a_1Ta^k+a_2T^k+a_3m+a_4RH+a_5I+a_0W

x=zip(oa['oa'][1:],zt[z_no][:-1],-1.*m[z_no][1:],rh[z_no][1:],q[z_no][1:])
y=zt[z_no][1:]

regr.fit(x, y)

lin.fit(x, y)

x2=zip(oa['oa'][1:],zt[z_no][:-1],m[z_no][1:],rh[z_no][1:],q[z_no][1:],[1]*(len(oa)-1))
res = lsq_linear(x2,y,bounds=([0, 0,  -np.inf, 0,0,-np.inf], [np.inf,0.99, 0,np.inf,np.inf,np.inf]))
cs=res['x']

print regr.coef_

print lin.coef_

print cs

a1=cs[0]

a2=cs[1]

a3=cs[2]

a4=cs[3]

a5=cs[4]

a0=cs[5] 


c1=regr.coef_[0]

c2=regr.coef_[1]

c3=-regr.coef_[2]

c4=regr.coef_[3]

c5=regr.coef_[4]

c0=regr.intercept_ 

y=regr.predict(x)
y1=lin.predict(x)

print len(y)

t_pre1=zt[z_no].iloc[0]
t_pre2=zt[z_no].iloc[0]
t_pre3=zt[z_no].iloc[0]
t_pre4=zt[z_no].iloc[0]

y_pred1=[t_pre1]
y_pred2=[t_pre2]
y_pred3=[t_pre3]
y_pred4=[t_pre4]

for i in range(1,len(rh)):
    t_pre1=a1*oa['oa'].iloc[i]+a2*t_pre1+a3*m[z_no].iloc[i]+a4*rh[z_no].iloc[i]+a5*q[z_no].iloc[i]+a0
    t_pre2=c1*oa['oa'].iloc[i]+c2*t_pre2+c3*m[z_no].iloc[i]+c4*rh[z_no].iloc[i]+c5*q[z_no].iloc[i]+c0
    t_pre3=a1*oa['oa'].iloc[i]+a2*zt[z_no].iloc[i-1]+a3*m[z_no].iloc[i]+a4*rh[z_no].iloc[i]+a5*q[z_no].iloc[i]+a0
    t_pre4=c1*oa['oa'].iloc[i]+c2*zt[z_no].iloc[i-1]+c3*m[z_no].iloc[i]+c4*rh[z_no].iloc[i]+c5*q[z_no].iloc[i]+c0
    y_pred1.append(t_pre1)
    y_pred2.append(t_pre2)
    y_pred3.append(t_pre3)
    y_pred4.append(t_pre4)


xx=np.arange(len(y_pred1))

plt.plot(xx,y_pred1, label='prediction1')
#plt.plot(xx[1:],y, label='prediction2')
#plt.plot(xx,y_pred2, label='prediction2')
#plt.plot(xx[1:],y1, label='prediction2 fb')
plt.plot(xx,zt[z_no], label='real')
plt.legend()
#plt.xticks(xtic,xlab,rotation=45)
plt.xlim(xx[0],xx[-1])
plt.savefig('result/'+z_no+'.png',bbox_inches = 'tight',pad_inches = 0.1)
plt.show()

coefs=collections.OrderedDict()
coefs['a0']=a0
coefs['a1']=a1
coefs['a2']=a2
coefs['a3']=a3
coefs['a4']=a4
coefs['a5']=a5

with open('result/'+z_no+'.config','w') as outfile:
      json.dump(coefs,outfile, sort_keys=True, indent=2) 


