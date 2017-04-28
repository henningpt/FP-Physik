import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import scipy.constants as sc
import uncertainties.unumpy as unp
# daten laden
# 10^-11 , 0.1

I1, T1 = np.genfromtxt("v48m1.txt", unpack=True)
# 10^-11 , 0.3
I2, T2 = np.genfromtxt("v48m12.txt", unpack=True)
# 10^-11 , 1.0
I3, T3 = np.genfromtxt("v48m13.txt", unpack=True)
# 10^-10, 0.3
I4, T4 = np.genfromtxt("v48m14.txt", unpack=True)

I5, T5 = np.genfromtxt("v48m15.txt", unpack=True)
# groessen definieren
I1 *= 0.1 * 10**(-11)
I2 *= 0.3 * 10**(-11)
I3 *= 1.0 * 10**(-11)
I4 *= 0.3 * 10**(-10)
I5 *= 0.3 * 10**(-10)

# haenge daten aneinander
I12 = np.append(I1, I2)
I34 = np.append(I3, I4)
I14 = np.append(I12, I34)
I = np.append(I14, I5)

T12 = np.append(T1, T2)
T34 = np.append(T3, T4)
T14 = np.append(T12, T34)
T = np.append(T14, T5)

# funktionen


def lin(Tinv, intercept, slope):
    return(intercept + Tinv * slope)


# rechnung

Ic = I[:10]
Tc = T[:10]
params, cov = curve_fit(lin, 1/Tc, np.log(Ic))
errors = np.sqrt(np.diag(cov))

interc = unp.uarray(params[0], errors[0])
slope  = unp.uarray(params[1], errors[1])
# ground = lin(T, *params)


# plotten
x_plot = np.linspace(-61.8, -35.0, 1000)
plt.plot(1/Tc, np.log(Ic), 'rx')
# plt.plot(1/T, np.log(I) - ground, 'gx')
plt.plot(1/x_plot, lin(1/x_plot, *params), 'b')
plt.show()


# ausgabe

# W mit methode 1
print("\nslope: ", slope)
print("\nintercept: ", interc)
