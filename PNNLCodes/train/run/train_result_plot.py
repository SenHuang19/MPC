import pandas as pd
import matplotlib.pyplot as plt
import numpy as np






df1 = pd.read_csv('r2.csv', header=None)

df2 = pd.read_csv('me.csv', header=None)

df3 = pd.read_csv('r2_idea.csv', header=None)

df4 = pd.read_csv('me_idea.csv', header=None)

df5 = pd.read_csv('a2.csv', header=None)



df1=df1.sort_values(by=[0])


df2=df2.sort_values(by=[0])


df3=df3.sort_values(by=[0])


df4=df4.sort_values(by=[0])

df5=df5.sort_values(by=[0])



for ii in range(16):

     xx=np.arange(len(df1[0]))
     fig, ax1 = plt.subplots(3,1)


     ax1[0].plot(df1[0], df1[ii+1], label='long-term',marker='o',color='r')
     ax1[0].plot(df1[0], df3[ii+1], label='short-term',marker='x',color='b')
     
     ax1[1].plot(df2[0], df2[ii+1], label='long-term',marker='o',color='r')
     ax1[1].plot(df2[0], df4[ii+1], label='short-term',marker='x',color='b')

     ax1[2].plot(df2[0], df5[ii+1],marker='.')

     ax1[1].set_xlabel('Discrete Interval [min]')
     # ax1.set_ylim(df3[ii+1].min()-0.01,df1[ii+1].max()+0.01)
     # ax2.set_ylim(df4[ii+1].min()-0.01,df2[ii+1].max()+0.01)
     ax1[0].set_ylabel('$R^2$')
     ax1[1].set_ylabel('$RMSE$ [$^oC$]')
     ax1[2].set_ylabel('$a2$')
     ax1[0].legend()

     plt.savefig(str(ii)+'_result.png',bbox_inches = 'tight',pad_inches = 0.1)





