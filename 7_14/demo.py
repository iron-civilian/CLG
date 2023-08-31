import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import pickle
from scipy.optimize import curve_fit as cf
import sys
#%matplotlib notebook
import time


T_arr=np.linspace(0.8,2.94285714,26)

ext_arr=list(range(1,27))
data_arr=[]
for i in ext_arr:
    data_arr=data_arr+list(np.loadtxt(f'data_{i}.csv'))
data_arr=np.array(data_arr)
data_arr=data_arr.reshape((-1,3))

m_arr=data_arr[:,0]
m2_arr=data_arr[:,1]
m4_arr=data_arr[:,2]

np.savetxt('7_14.txt',data_arr)
#plt.plot(T_arr,m_arr,marker='o')
#plt.show()


