import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.constants as sc
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

# funktionen


def logj(Tinv, a, W):
    return(a - Tinv * W/(sc.k))


# rechnung

Ic = I[:16]
Tc = T[:16]

Ilog = np.log(I)
Tinv = 1/T

# params, cov = curve_fit(j, Tc, Ic)


# plotten
# x_plot = np.linspace(-61.8, 63.8, 1000)
plt.plot(T, I, 'rx')

# plt.plot(x_plot, j(x_plot,*params), 'b')
plt.show()


# ausgabe
print(Tinv)
