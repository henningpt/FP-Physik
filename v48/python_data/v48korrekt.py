import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import scipy.constants as sc
import uncertainties.unumpy as unp
from uncertainties import ufloat
import scipy.integrate as integrate

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
T = T + 273.15

b = abs(T[0] - T[len(T) - 1]) / (len(T) - 1) / 60
# funktionen


def kreisf(r):
    return(np.pi * r**2)


def nonegatives(arr):
    l = len(arr)
    for i in range(0, l):
        if (arr[i] < 0):
            arr[i] = I[0]
    return(arr)


def lin(Tinv, intercept, slope):
    return(intercept + Tinv * slope)


def tau0(Tmax, W, b):
    return(sc.k * Tmax**2 * np.e**(-W / (sc.k * Tmax)) / (W * b))


def efct(T, a, b):
    return(a * np.e**(T * b))


# rechnung
A = kreisf(2.5*10**(-3))

# untergrundfit + abziehen
Iu = np.append(I[9:11], I[28:32])
Tu = np.append(T[9:11], T[28:32])

params_u, covariance_u = curve_fit(lin, Tu, Iu)
eparams_u = unp.uarray(params_u, np.sqrt(np.diag(covariance_u)))
slope_u = eparams_u[1]
interc_u = eparams_u[0]

Ik = I - lin(T, *params_u)
Ik2 = Ik[I > Ik]
Tk2 = T[I > Ik]
Ir = Ik2[Ik2 >= 0]
Tr = Tk2[Ik2 >= 0]

# Anfangsfit
lin_fits = 3
lin_fite = 10
Ic = Ir[lin_fits:lin_fite]
jc = Ic/A  # stromdichte
Tc = Tr[lin_fits:lin_fite]
params_lin, cov_lin = curve_fit(lin, 1/Tc, np.log(Ic))
errors_lin = np.sqrt(np.diag(cov_lin))
interc_lin = ufloat(params_lin[0], errors_lin[0])
slope_lin = ufloat(params_lin[1], errors_lin[1])
W = -slope_lin * sc.k / sc.elementary_charge


# plotten
# werte auswahl etc.
plt.figure(1)  # T, I 1
plt.plot(T, I, 'rx', label=r'$\mathrm{unkorr.} \ \mathrm{Werte}$')
plt.plot(T, lin(T, *params_u), label=r'$\mathrm{Untergrund-Fit}$')
plt.plot(T, Ik, 'gx',  label=r'$\mathrm{korr.} \ \mathrm{Werte}$')
plt.plot(Tr, Ir, 'x', color=(1.0, 0.0, 1.0),
         label=r'$\mathrm{pos.korr.} \ \mathrm{Werte}$')

# plt.plot(Tc, Ic, 'gx', label=r'$\mathrm{Werte} \ \mathrm{fuer} \ \mathrm{Anfangskurve}$')
plt.xlabel(r'$T / \mathrm{K}$')
plt.ylabel(r'$I / \mathrm{A}$')
plt.legend(loc='best')

# anfangsplot
plt.figure(2)
x_plot = np.linspace(np.min(Tc), np.max(Tc), 1000)
plt.plot(1/Tc, np.log(Ic), 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
plt.plot(1/x_plot, lin(1/x_plot, *params_lin), 'b', label=r'$\mathrm{linearer} \ \mathrm{Fit}$')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln(j \cdot \mathrm{m^2} \ / \ \mathrm{A}) $')
plt.legend(loc='best')
plt.show()


# ausgabe
#lin fit ergebnis W
print("Aktivierungsenergie W: ", W)
