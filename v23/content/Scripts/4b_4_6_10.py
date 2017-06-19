import numpy as np
import math
# fit arbitrary functions
from scipy.optimize import curve_fit

A=np.genfromtxt('4b_4_6_10.dat', unpack=True)
X=[]

#Parameter um peaks zu pruefen, Anzahl der Werte die kleiner sind als A_1_i+5
m=0

#Filter Bauartbedingte Resonanzen der Soundkarte
f_1=310
f_2=2130
f=800
f_delta=20

#minimale Peak hoehe
f_abs=0.3
f_abs_3000=2
#abstand zwischen zwei peaks
r=5

#quotientenparameter
q=1

for i in range(len(A[1,:])-r+1):
    for j in range(r):
        if A[1,i+r/2-1] > q*A[1,i+j]:
            m = m+1 
    if A[0,i+r/2-1]>3000:
        if m>= r-1  and A[0,i+r/2-1]>f_1+f_delta and A[1,i+r/2-1]>f_abs:
            X = X+[A[:,i+r/2-1]]
        elif m>= r-1 and A[0,i+r/2-1]<f_1-f_delta and A[1,i+r/2-1]>f_abs:
            X = X+[A[:,i+r/2-1]]                
    m=0

Y=np.column_stack(X)
print(Y)
