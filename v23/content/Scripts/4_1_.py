import numpy as np 
import math

A=np.genfromtxt('4_1_peaks.txt', unpack=True)
Out=A
for i in range(7):
    for j in range(13):
        Out[i,j]=A[i,j+1]-A[i,j]
np.savetxt('4_1_diff.txt', np.column_stack(Out),delimiter='&',fmt='%3d',newline='\ ')
mean0=Out[0,0:13]
mean1=Out[1,0:10]
mean2=Out[2,0:9]
mean3=Out[3,0:8]
mean4=Out[4,0:6]
mean5=Out[5,0:4]
mean6=Out[6,0:4]
L=[0.600,0.525,0.450,0.375,0.300,0.225,0.150]
mean=[np.mean(mean0),np.mean(mean1),np.mean(mean2),np.mean(mean3),np.mean(mean4),np.mean(mean5),np.mean(mean6)]
std=[np.std(mean0),np.std(mean1),np.std(mean2),np.std(mean3),np.std(mean4),np.std(mean5),np.std(mean6)]

c=[0,0,0,0,0,0,0]
for i in range(7):
    c[i]=2*L[i]*mean[i] 

c_err=[0,0,0,0,0,0,0]
for i in range(7):
    c_err[i]=2*L[i]*std[i] 

a=0
for i in range(7):
    a=a+c_err[i]*c_err[i]

c_std=math.sqrt(a)/7
c_mean=np.mean(c)
print(c)
print(c_err)
print(c_mean,'+-',c_std)