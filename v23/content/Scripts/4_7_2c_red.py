import numpy as np
import math
# fit arbitrary functions
from scipy.optimize import curve_fit

A=np.genfromtxt('4b_7_2c.dat', unpack=True)
X=[]
print(len(A[0,:]))
print(A[1,0])
print(A[1,1])
#Parameter um peaks zu pruefen, Anzahl der Werte die kleiner sind als A_1_i+5
m=0

#Filter Bauartbedingte Resonanzen der Soundkarte
f_1=820
f_2=2130
f=800
f_delta=20

#minimale Peak hoehe
f_abs=0.2
f_abs_3000=0
#abstand zwischen zwei peaks
r=9

#quotientenparameter
q=1

for i in range(len(A[1,:])-r+1):
    for j in range(r):
        if A[1,i+r/2-1] > q*A[1,i+j]:
            m = m+1
    #Ausnahmen   
    #if A[0,i+r/2]==3380:
    #    X=X+[A[:,i+r/2-1]]
    #if A[0,i+r/2]==10690:
    #   X=X+[A[:,i+r/2-1]]    
    #if A[0,i+r/2]==10970:
    #  X=X+[A[:,i+r/2-1]]   
    #if A[0,i+r/2]==10040:
    #   X=X+[A[:,i+r/2-1]]          
    #regulaerer Ablauf
    if A[0,i+r/2-1]<8500 and not A[0,i+r/2-1]==800 and not A[0,i+r/2-1]==150 and not A[0,i+r/2-1]==280 :
        if m>= r-1  and A[0,i+r/2-1]>f_1+f_delta and A[1,i+r/2-1]>f_abs:
            X = X+[A[:,i+r/2-1]]
        elif m>= r-1 and A[0,i+r/2-1]<f_1-f_delta and A[1,i+r/2-1]>f_abs:
            X = X+[A[:,i+r/2-1]]      
    m=0

Y=np.column_stack(X)

k=np.arange(0,len(Y[0,:])+10)

k=k* math.pi /0.6
#print('k',k[1:3])
#k.itemset(3,111)
#print('k ',k)
#print(len(k))
#print('type k ',type(k))
#print('type Y ',type(Y))
print('Y ',Y)
#print('k.T ',k.T)
#print('Y.T ',Y[0,:].T)




# prepare plot

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 15)
plt.rcParams['font.size'] = 16

test=np.zeros(41)


#plt.plot(k[1:len(Y[0,:])+1],Y[0,:],'bx')

plt.plot(k[1:11],Y[0,0:10],'rx')
plt.plot(k[2:11],Y[0,10:19][::-1],'rx')
plt.plot(k[1:10],Y[0,19:28],'rx')

plt.xlabel('Wellenzahl k in 1/m')
plt.ylabel('Resonanzfrequenzen in Hz')
plt.legend(loc="best")
#plt.show()
plt.savefig('4b_7_2c.jpg')
#np.savetxt('4_2_peaks.txt', np.column_stack(X).T,delimiter='&',fmt='%3d',newline='\ xx ')

#print('c =', params[0]*2*math.pi, '+-', errors[0]*2*math.pi)
