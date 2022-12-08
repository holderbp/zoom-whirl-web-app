import math
from math import *
import numpy as np
from numpy import *
from numpy.linalg import *
import random
import matplotlib.pyplot as plt
J=1 #interaction energy
N=81 #number of cells
D=int(sqrt(N)) #dimension of lattice
nmc=1000000 #number of monte carlo steps for each temperature value
nmceq=100000 #number of states that are counted for expectation value
#k=1.380649*(10**-23) #boltzmann constant in m^2 kg s^-2 K^-1
k=1
Tlist=[]
Mlist=[]

def getmagnetization(x):
    M=0
    for i in range (0, N):
        M=M+latt[i]
    return M

def randspin():
    a=random.random()
    if a < 0.5:
        return -1
    else:
        return 1
    
lattup=np.zeros(N)
lattup=lattup.astype(int)
A=np.random.rand(N)
B=(A<0.5)
C=B.astype(int)

def getleftneighbor(cell):
    index = (cell - 1)%N
    return index
def getrightneighbor(cell):
    index = (cell + 1)%N
    return index
def getbottomneighbor(cell):
    index = (cell - D)%N
    return index
def gettopneighbor(cell):
    index = (cell + D)%N
    return index
def getleftneighborvalue(cell):
    value=latt[getleftneighbor(cell)]
    return value
def getrightneighborvalue(cell):
    value=latt[getrightneighbor(cell)]
    return value
def gettopneighborvalue(cell):
    value=latt[gettopneighbor(cell)]
    return value
def getbottomneighborvalue(cell):
    value=latt[getbottomneighbor(cell)]
    return value
def flipup(x):
    x=lattup
    return x

    
#def getaccprob(x):
    #possible dE values: -8*J, -4*J, 0, 4*J, 8*J
    #possible accprob values: math.exp(-4*beta*J)), math.exp(-8*beta*J))
Tarray=[]
for i in range (15, 50):
    Tarray.append(i/10)
for T in Tarray:
    N0=0
    N1=0
    N2=0
    latt=2*C-1
    print(latt)
    if T==Tarray[0]:
        print('first magnetization is: ', getmagnetization(latt))
    beta=1/(k*T)
    acc1=math.exp(-4*beta*J)
    acc2=math.exp(-8*beta*J)
    M=0.0
    Madd=0.0
    for i in range (0, N):
        M=M+latt[i]
    for j in range (0, nmc):
        celltoflip=random.randint(0, N-1)
        oldspin=latt[celltoflip]
        dE=2.0*J*latt[celltoflip]*(latt[gettopneighbor(celltoflip)]+latt[getbottomneighbor(celltoflip)]+latt[getrightneighbor(celltoflip)]+latt[getleftneighbor(celltoflip)])
        if dE<=0:
            latt[celltoflip] = -latt[celltoflip]
            N0=N0+1
        else:
            #if dE>3.9*J and dE<4.1J:
            if dE==4.0*J:
                rand=random.random()
                if rand < acc1:
                    latt[celltoflip] = -latt[celltoflip]
                    N1=N1+1
            #if dE<8.1*J and dE>7.9*J:
            if dE==8.0*J:
                rand=random.random()
                if rand<acc2:
                    latt[celltoflip] = -latt[celltoflip]
                    N2=N2+1
        if oldspin != latt[celltoflip]:
            M=M+2.0*latt[celltoflip]
        if j>nmc-nmceq:
            Madd=Madd+M
    Madd=Madd/nmceq/D**2
    Mlist.append(abs(Madd))
    Tlist.append(T)
    T=T+0.1
    print(latt)
    print('new magnetization is:', getmagnetization(latt))
    print(N0, N1, N2)
plt.plot(Tlist, Mlist, 'o')
plt.show()
#print('new energy is:', E)
