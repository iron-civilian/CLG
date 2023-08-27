import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import pickle
from scipy.optimize import curve_fit as cf
import sys
#%matplotlib notebook
import time
 

kB=1
J=1
L=8#int(sys.argv[1])
T=2#int(sys.argv[2])
e=6#int(sys.argv[2])
T1=0.2#float(sys.argv[3])
T2=3#float(sys.argv[4])
T_N=6#int(sys.argv[5])
Lx=L
Ly=2*L
N_sites=Lx*Ly
#N_steps=10**e

start = time.time()

def nbr2D(L):
    N_sites=Lx*Ly
    nbrarr=np.zeros((N_sites,4),dtype=int)
    
    for i in range(N_sites):
        k1=i//Ly
        k2=i%Ly
        if 1<=k1<=Lx-2:
            nbrarr[i][2]=(k1-1)*Ly+k2 #up
            nbrarr[i][3]=(k1+1)*Ly+k2 #down
        else:
            a=(k1==0)
            nbrarr[i][2]=Ly*(Lx-2+a)+k2 #up
            nbrarr[i][3]=Ly*a+k2 #down
        if 1<=k2<=L-2:
            nbrarr[i][0]=i-1 #left
            nbrarr[i][1]=i+1 #right
        else:
            b=(k2==0)
            nbrarr[i][0]=(i-1+Ly*b)
            nbrarr[i][1]=(i+1-Ly*(1-b))
    return nbrarr
nbrarr=nbr2D(L)

'''
@jit(nopython=True)
def get_ensemble_avg(T,N_steps,N_ensemble,L):
    Sum_m=np.zeros(N_steps)
    Sum_E=np.zeros(N_steps)
    
    for i in range(N_ensemble):
        m_arr,E_arr=MC_update(T,N_steps,L)
        Sum_m=Sum_m+m_arr
        Sum_E=Sum_E+E_arr
    return Sum_m/N_ensemble,Sum_E/N_ensemble
'''

@jit(nopython=True)
def MC_update(state_arr,E,T): 
    beta=1/T

    #state_arr=np.ones(N_sites)
    #beta=1/T #k_B=1
    i1=np.random.randint(N_sites) #site chosen
    i2=np.random.randint(4) #neighbour chosen
    r=np.random.random()
    delE=0
    for k in range(4):
    	delE=delE+state_arr[nbrarr[i1][k]]
    for k in range(4):
        delE=delE-state_arr[nbrarr[nbrarr[i1][i2]][k]]
    
    delE=J*delE
    
    if (state_arr[i1]==1 and state_arr[nbrarr[i1][i2]]==0) and delE<0 or r<np.exp(-beta*delE):
        #m=m-2*state_arr[i]/N_sites
        E=E+delE
        state_arr[i1]=0
        state_arr[nbrarr[i1][i2]]=1
    return E

'''    
@jit(nopython=True)
def eq_fit(t,tau):
    return np.exp(-t/tau)
'''


T_arr=np.linspace(T1,T2,T_N)

'''
tau_arr=np.zeros(T_N)
get_ensemble_avg(2,10,10,4)
N_range=np.arange(1,2501,1)
for i in range(T_N):
    m_arr,E_arr=get_ensemble_avg(T_arr[i],2500,1000,L)
    popt,pcov=cf(eq_fit,N_range,m_arr)
    tau_arr[i]=popt[0]
'''


'''
m_mean_arr=np.zeros(T_N)
m2_mean_arr=np.zeros(T_N)
m4_mean_arr=np.zeros(T_N)
for i in range(T_N):
	m_arr,E_arr=get_ensemble_avg(T_arr[i],int(tau_arr[i])*10+10**e,1,L)
	m_mean_arr[i]=np.mean(m_arr[-10**e:])
	m2_mean_arr[i]=np.mean(m_arr[-10**e:]**2)
	m4_mean_arr[i]=np.mean(m_arr[-10**e:]**4)
'''
#@jit(nopython=True)
def get_data(T,N_ensemble):
    #datafile=open('data.txt','a')
    state_arr=np.array((L*[1]+L*[0])*L)
    #np.savetxt(datafile,state_arr)

    Trlax=10000

    
    E=-J*((L*L)+L*(L-1))
    E_mean=0
    E2_mean=0
    E4_mean=0
    for i in range(Trlax*N_sites):
        E=MC_update(state_arr,E,T)
        #E4_mean+=E*E*E*E
    for i in range(N_ensemble*N_sites):
        E=MC_update(state_arr,E,T)
        E_mean+=E
        E2_mean+=E*E
        E4_mean+=E*E*E*E
    return E_mean/N_ensemble,E2_mean/N_ensemble,E4_mean/N_ensemble


'''
for i in range(T_N):
    m,m2,m4=get_data(T_arr[i])
    m_mean_arr[i]=m
    m2_mean_arr[i]=m2
    m4_mean_arr[i]=m4
'''
E_mean_arr=np.zeros(T_N)
E2_mean_arr=np.zeros(T_N)
E4_mean_arr=np.zeros(T_N)

for i in range(T_N):
    E_,E2_,E4_=get_data(T_arr[i],10**e)
    E_mean_arr[i]=E_
    E2_mean_arr[i]=E2_
    E4_mean_arr[i]=E4_
    

plt.plot(T_arr,E_mean_arr)
plt.show()

end = time.time()

print(f'runtime = {end-start} sec')
