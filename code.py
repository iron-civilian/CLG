import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import pickle
from scipy.optimize import curve_fit as cf
import sys
#%matplotlib notebook
import time
 

kB=1
J=2
L=int(sys.argv[1])
T=float(sys.argv[2])
e=int(sys.argv[3])

Lx=L
Ly=2*L
N_sites=Lx*Ly
N_steps=10**e

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
        if 1<=k2<=Ly-2:
            nbrarr[i][0]=i-1 #left
            nbrarr[i][1]=i+1 #right
        else:
            b=(k2==0)
            nbrarr[i][0]=(i-1+Ly*b)
            nbrarr[i][1]=(i+1-Ly*(1-b))
    return nbrarr
nbrarr=nbr2D(L)

@jit(nopython=True)
def ord_param(state_arr):
    m=0
    for j in range(Ly):
        rho_y=0
        for i in range(Lx):
            rho_y+=state_arr[i*Ly+j]/Lx
        m+=np.abs(rho_y-0.5)
    return 2*m/Ly


@jit(nopython=True)
def MC_update(state_arr,E,T): 
    beta=1/T

    i1=np.random.randint(N_sites) #site chosen
    
    if state_arr[i1]==1:
        i2=np.random.randint(4)
        if state_arr[nbrarr[i1][i2]]==0:
            delE=0
            for k in range(4):
                delE=delE+state_arr[nbrarr[i1][k]]
            for k in range(4):
                delE=delE-state_arr[nbrarr[nbrarr[i1][i2]][k]]
            delE=2*J*(delE+1)
            if (delE<=0 or np.random.random()<np.exp(-beta*delE)):
                E+=delE
                state_arr[i1]=0
                state_arr[nbrarr[i1][i2]]=1
    return E  

@jit(nopython=True)
def H(state_arr):
    eng=0
    for i in range(N_sites):
        for k in range(4):
            eng+=state_arr[i]*state_arr[nbrarr[i][k]]  
    return -J*eng  

@jit(nopython=True)
def get_data(T):
    #datafile=open('data.txt','a')
    
    pos_arr=np.arange(0,N_sites,1)
    np.random.shuffle(pos_arr)
    pos_arr=pos_arr[:N_sites//2]
    state_arr=np.zeros(N_sites)
    for i in range(N_sites//2):
        state_arr[pos_arr[i]]=1
    
    #state_arr=np.array((L*[1]+L*[0])*L)
    #np.savetxt(datafile,state_arr)

    Trlax=10000000
    
    E=H(state_arr)
    #E_mean=0
    #E2_mean=0
    #E4_mean=0
    m_mean=0
    m2_mean=0
    m4_mean=0
    
    for i in range(Trlax*N_sites):
        E=MC_update(state_arr,E,T)
    
    #return state_arr
    
    for i in range(N_steps):
        for j in range(N_sites):
            E=MC_update(state_arr,E,T)
        m=ord_param(state_arr)
        m_mean+=m
        m2_mean+=m*m
        m4_mean+=m*m*m*m
 
    return m_mean/(N_steps),m2_mean/(N_steps),m4_mean/(N_steps)
    

data_T=get_data(T)
end = time.time()


#plt.imshow(state_arr.reshape(Lx,Ly))
#plt.show()

#data_T=get_data(T)
np.savetxt(f'data_{T}.csv',data_T,delimiter=',')




print(f'runtime = {end-start} sec')
