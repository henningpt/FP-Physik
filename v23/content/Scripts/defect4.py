import numpy as np
import math
# fit arbitrary functions
from scipy.optimize import curve_fit

A=np.genfromtxt('4b_7_12x50.dat', unpack=True)
B=np.genfromtxt('4b_7_2a.dat', unpack=True)

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (20,10)
plt.rcParams['font.size'] = 16

plt.plot(A[0,:],A[1,:],'b-',linewidth=0.5,label='ohne Defekt')

plt.plot(B[0,:],B[1,:],'r-',linewidth=0.5,label='mit Defekt')

plt.ylabel('Amplitude')
plt.xlabel('Frequenz in Hz')
plt.legend(loc="best")

plt.savefig('defect4.jpg')