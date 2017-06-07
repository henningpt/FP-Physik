import numpy as np

#Messwerte einlesen
Y=np.genfromtxt('v18_3.txt' , unpack=True)

X=np.arange(0,len(Y))

#Gaussfit
import array
from scipy.optimize import curve_fit

A=np.linspace(-10,10,21)
print(A)
print(A**2)

def f(x,a,b):
    return np.exp(-a*(x-b)**2)

#Plot Vorbereiten
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

plt.plot(X,Y[:],'rx')
plt.show