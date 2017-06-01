
import numpy as np
from scipy.optimize import curve_fit

WERTE=np.genfromtxt('v18_1.txt' , unpack=True)

Y=np.array([825,1382,1491,2611])

X=np.array([295,411,443,678])

A=np.array([121,244,295,344,411,443,678,688,778,867,964,1005,1085,1112,1299,1408,1457])
print(A)
def f(x,a,b):
    return a*x+b

params,covariance=curve_fit(f,X,Y)
print(f(A,*params))

#Funktionen
#Gauss
def gauss(a,s):
 X2=np.linspace(a-s,a+s,2*s+1)
# print('X2',X2)
 Y2=WERTE[a-s:a+s+1]
 #print('Y2',Y2)
 def g(x,p1,p2,p3,p4):
    return p4*np.e**(-p1*(x-p2)**2) + p3

 params2,covariance2 = curve_fit(g,X2,Y2, p0=[1,a,0,1])
 
 #Plot Vorbereiten
 import matplotlib.pyplot as plt
 plt.rcParams['figure.figsize'] = (10, 8)
 plt.rcParams['font.size'] = 16
 
 plt.plot(X2,Y2,'rx')
 plt.plot(X2,g(X2,*params2))
 #plt.show() 
 
 return params2

print('Parameter',gauss(825,20))

