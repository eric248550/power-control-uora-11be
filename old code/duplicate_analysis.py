import numpy as np
import matplotlib.pyplot as plt
import math
import random

# for the access success probability, Ps

# initial values
R=15 # the total RUs in each slot
K=[1,5,10] # the number of reservation after successful contention
Lmax=5 # the retry limit
TXOP=5.681 # transmission opportunity
T_TF=100*(10**(-6)) # period of trigger frame
T_SIFS=16*(10**(-6)) # period of SIFS
T_PIFS=25*(10**(-6)) # period of PIFS
T_RARU=5484*(10**(-6)) # period of RA-RU
T_MBA=40*(10**(-6)) # period of MBA
delta=[0,2,5] # the period of time for reserving one slot
P_PLR=[0,0.1,0.3,0.5] # the packet loss rate

def OCW(n):
    if n==1:
        return 7.0
    elif n==2:
        return 15.0
    else:
        return 31.0

iii=0

def analytical_model(M,K,delta=0):
    Imax=200 # initialize Imax

    Mi=[0 for _ in range(0,Imax+1)]
    Ri=[0 for _ in range(0,Imax+1)]
    M_iS=[0 for _ in range(0,Imax+1)]
    M_iF=[0 for _ in range(0,Imax+1)]

    M_in=[[0]*(Lmax+1) for _ in range(0,Imax+1)]
    M_iSn=[[0]*(Lmax+1) for _ in range(0,Imax+1)]
    M_iFn=[[0]*(Lmax+1) for _ in range(0,Imax+1)]
    
    for i in range(1,Imax+1):
        # get Ri
        if i==1:
            Ri[i]=R
        elif i>1 and i<=K:
            Ri[i]=Ri[i-1]-M_iS[i-1]
        elif i>=(K+1):
            Ri[i]=Ri[i-1]-M_iS[i-1]+M_iS[i-K]

        for n in range(1,Lmax+1):
            # initial condition
            a_j_i=0
            if n==1:
                if i==1 and Ri[i]>=OCW(n):
                    a_j_i=1
                elif i==1 and Ri[i]<OCW(n):
                    a_j_i=(float)((Ri[i]+1)/(OCW(n)+1))
                elif i>1:
                    temp=0
                    for c in range(1,i):
                        temp=temp+Ri[c]
                    if (OCW(n)-temp)>=Ri[i]:
                        a_j_i=(float)((Ri[i])/(OCW(n)+1))
                    elif (OCW(n)-temp)<=0:
                        a_j_i=0
                    else:
                        a_j_i=(float)((OCW(n)-temp)/(OCW(n)+1))
                M_in[i][n]=a_j_i*M
            elif n<=i:
                for k in range(1,i):
                    if k==i-1 and Ri[i]<OCW(n):
                        a_j_i=(float)((Ri[i]+1)/(OCW(n)+1))
                    elif k==i-1 and Ri[i]>=OCW(n):
                        a_j_i=1
                    elif k<i-1:
                        temp=0
                        for c in range(k+1,i):
                            temp=temp+Ri[c]
                        if OCW(n)-temp>=Ri[i]:
                            a_j_i=(float)((Ri[i])/(OCW(n)+1))
                        elif (OCW(n)-temp)<=0:
                            a_j_i=0
                        else:
                            a_j_i=(float)((OCW(n)-temp)/(OCW(n)+1))
                    M_in[i][n]=M_in[i][n]+a_j_i*M_iFn[k][n-1]
            Mi[i]=Mi[i]+M_in[i][n]
        # statistics
        for n in range(1,Lmax+1):
            if Ri[i]!=0:
                M_iSn[i][n]=float(M_in[i][n]*np.exp(-Mi[i]/Ri[i]))
            M_iFn[i][n]=float(M_in[i][n]-M_iSn[i][n])
            M_iS[i]=M_iS[i]+M_iSn[i][n]

    # get actual Imax
    Nn_index=0
    Imax=0
    for n in range(1,Lmax+1):
        Sum=0.0
        Nn=0
        for c in range(Nn_index+1,300+1):
            if Sum<=OCW(n):
                Nn=Nn+1
                Nn_index=Nn_index+1
                Sum=Sum+Ri[c]
            else:
                break
        Imax=Imax+Nn
    Imax=Imax+(K-1)

    # get Da
    total_success=0.0
    for i in range(1,Imax+1):
        total_success=total_success+M_iS[i]*K
    total_delay=0.0
    for i in range(1,Imax+1):
        temp=0.0
        for t in range(i,(i+K-1)+1):
            temp=temp+t
        total_delay=total_delay+M_iS[i]*temp

    Da=total_delay*TXOP/total_success
    return Da

def simulation(M,K):
    return 0
    samples=1000
    total_samples=0.0

    for s in range(samples):
        total_success=0.0
        total_packets=M*K
        
        # get Ps
        total_samples=total_samples+total_success/total_packets
    return total_samples/samples

if __name__=="__main__":
    # declare the array
    k1_ana=[]
    k1_sim=[]
    k5_ana=[]
    k5_sim=[]
    k10_ana=[]
    k10_sim=[]

    M_ana=np.arange(0, 101, 10) # the number of stations
    M_ana[0]=1
    for m in M_ana:
        print(m)
        k1_ana.append(analytical_model(m,K[0]))
        k5_ana.append(analytical_model(m,K[1]))
        k10_ana.append(analytical_model(m,K[2]))
    
    M_sim=np.arange(0, 101, 10) # the number of stations
    M_sim[0]=1
    for m in M_sim:
        k1_sim.append(simulation(m,K[0]))
        k5_sim.append(simulation(m,K[1]))
        k10_sim.append(simulation(m,K[2]))

    # plot the figure
    plt.plot(M_sim,k1_sim,linestyle='None',marker='+',markerfacecolor='white',color='blue',label='Sim K=1')
    plt.plot(M_ana,k1_ana,linestyle='dashdot',color='blue',label='Ana K=1')
    plt.plot(M_sim,k5_sim,linestyle='None',marker='^',markerfacecolor='white',color='red',label='Sim K=5')
    plt.plot(M_ana,k5_ana,linestyle='dashed',color='red',label='Ana K=5')
    plt.plot(M_sim,k10_sim,linestyle='None',marker='o',markerfacecolor='white',color='green',label='Sim K=10')
    plt.plot(M_ana,k10_ana,linestyle='solid',color='green',label='Ana K=10')

    plt.xlabel("number of stations (M)")
    plt.ylabel("access success probability")
    plt.xlim((0,100))
    plt.ylim((0,90))

    plt.legend()
    plt.show()