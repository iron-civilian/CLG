import numpy as np
import matplotlib.pyplot as plt

T_arr=np.linspace(0.8,2.94285714,26)
data_arr1=np.loadtxt('7_14.txt')
data_arr2=np.loadtxt('8_16.txt')
data_arr3=np.loadtxt('10_20.txt')

L_arr=[7,8,10]

m_arr1=data_arr1[:,0]
m_arr2=data_arr2[:,0]
m_arr3=data_arr3[:,0]


m2_arr1=data_arr1[:,1]
m2_arr2=data_arr2[:,1]
m2_arr3=data_arr3[:,1]


m4_arr1=data_arr1[:,2]
m4_arr2=data_arr2[:,2]
m4_arr3=data_arr3[:,2]

b1=m2_arr1**2/m4_arr1
b2=m2_arr2**2/m4_arr2
b3=m2_arr3**2/m4_arr3


plt.plot((T_arr-1.13)*7,m_arr1*7**(1/8),label='7_14')
plt.plot((T_arr-1.13)*8,m_arr2*8**(1/8),label='8_16')
plt.plot((T_arr-1.13)*10,m_arr3*10**(1/8),label='10_20')
plt.legend()
plt.show()




