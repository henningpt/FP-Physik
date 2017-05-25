import numpy as np
import math
# fit arbitrary functions
from scipy.optimize import curve_fit

A=np.genfromtxt('4_3_10mm.dat', unpack=True)
X=[]
print(len(A[0,:]))
print(A[1,0])
print(A[1,1])
#Parameter um peaks zu pruefen, Anzahl der Werte die kleiner sind als A_1_i+5
m=0

#minimale frequenz
f_min=400

#minimale Peak hoehe
f_abs=2

for i in range(len(A[1,:])-10):
    for j in range(11):
        if A[1,i+5] > A[1,i+5+j-5]:
            m = m+1
    if m>= 10  and A[0,i+5]>f_min and A[1,i+5]>f_abs:
        X = X+[A[:,i+5]]
    m=0

Y=np.column_stack(X)

k=np.arange(1,len(Y[0,:])+1)

k=k* math.pi /0.6
print(k)
# prepare plot

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

test=np.zeros(41)


#def f(x, a, b):
#    return a * x + b

# covariance is the covariance matrix

#errors = np.sqrt(np.diag(covariance))

#print('a =', params[0], '+-', errors[0])
#print('b =', params[1], '+-', errors[1])

#x_plot = np.linspace(1, k[len(Y[0,:])-1])

#plt.plot(x_plot, f(x_plot, *params), 'b-', label='f(k)=54.96 k + 20', linewidth=2)

plt.plot(k,np.column_stack(X)[0,:],'rx',label="Messdaten")
plt.xlabel('Wellenzahl k in 1/m')
plt.ylabel('Resonanzfrequenzen in Hz')
plt.legend(loc="best")
plt.show()
#np.savetxt('4_2_peaks.txt', np.column_stack(X).T,delimiter='&',fmt='%3d',newline='\ xx ')

#print('c =', params[0]*2*math.pi, '+-', errors[0]*2*math.pi)
