import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.constants as sc
import uncertainties.unumpy as unp

# daten laden
# 10^-11 , 3
I1, T1 = np.genfromtxt("2v48m1.txt", unpack=True)

# 10^-11 , 0.3
I2, T2 = np.genfromtxt("2v48m12.txt", unpack=True)

# 10^-11 , 1.0
I3, T3 = np.genfromtxt("2v48m13.txt", unpack=True)

# 10^-10, 0.3
I4, T4 = np.genfromtxt("2v48m14.txt", unpack=True)

# 10^-11 1.0
I5, T5 = np.genfromtxt("2v48m15.txt", unpack=True)

# groessen definieren
I1 *= 3 * 10**(-11)
I2 *= 0.3 * 10**(-11)
I3 *= 1.0 * 10**(-11)
I4 *= 0.3 * 10**(-10)
I5 *= 1.0 * 10**(-11)

# haenge daten aneinander
I12 = np.append(I1, I2)
I34 = np.append(I3, I4)
I14 = np.append(I12, I34)
I = np.append(I14, I5)

T12 = np.append(T1, T2)
T34 = np.append(T3, T4)
T14 = np.append(T12, T34)
T = np.append(T14, T5)
T = T + 273.15
# sortiere arrays
Isort = np.zeros(len(I))
Thelp = np.argsort(T)
for i in range(0, len(I)):
    Isort[i] = I[Thelp[i]]

Tsort = np.sort(T)


# funktionen


def lin(Tinv, intercept, slope):
    return(intercept - Tinv * slope)


# rechnung

Tinv = 1/Tsort
#Ic = np.delete(Isort[:20], 5)
#Tc = np.delete(Tsort[:20], 5)

Ic = Isort[:20]
Tc = Tsort[:20]

params, cov = curve_fit(lin, 1/Tc, np.log(Ic))
errors = np.sqrt(np.diag(cov))

interc = unp.uarray(params[0], errors[0])
slope = unp.uarray(params[1], errors[1])
W = slope * -sc.k / sc.elementary_charge


# plotten
x_plot = np.linspace(-61.8 + 273.15, -30.0 + 273.15, 1000)

plt.figure(1)
plt.xlabel('T / K')
plt.ylabel('I / A')
plt.plot(T, I, 'rx')
plt.savefig('ti2.pdf')
print(Ic)

plt.figure(2)
plt.plot(Tc, Ic, 'gx')
#plt.xlabel('')
#plt.ylabel('')
plt.plot(1/x_plot, lin(1/x_plot, *params), 'b')
plt.savefig('lin2.pdf')


# ausgabe
print("W: ", W)
