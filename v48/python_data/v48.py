48import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import scipy.constants as sc
import uncertainties.unumpy as unp
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
    return(a * np.e**(T * b

# rechnung
A = kreisf(2.5*10**(-3))


Iu = I[28:33]
Tu = I[28:33]

uparams, ucov = curve_fit(lin, Tu, Iu)
uerrors = np.sqrt(np.diag(ucov))
uinterc = unp.uarray(uparams[0], uerrors[0])
uslope = unp.uarray(uparams[1], uerrors[1])

Ik = np.zeros(len(I))
Ik[0:32] = I[0:32] - lin(T[0:32], *uparams)
Ik[32:] = I[32:]
# Ik = nonegatives(Ik)

Ic = I[:16]
jc = Ic/A  # stromdichte
Tc = T[:16]
params, cov = curve_fit(lin, 1/Tc, np.log(jc))
errors = np.sqrt(np.diag(cov))

interc = unp.uarray(params[0], errors[0])
slope = unp.uarray(params[1], errors[1])


W = slope * -sc.k / sc.elementary_charge

maxgrenz = 29
mingrenz = 4

print("Ik", Ik)
print("T,", T)
integ = integrate.simps(T[mingrenz:maxgrenz], Ik[mingrenz:maxgrenz])  # integration ueber strom
l = integ/Ik  # / b
print("integral: ", integ)
params2, cov2 = curve_fit(lin, 1/T[mingrenz:maxgrenz], np.log(l[mingrenz:maxgrenz]))
errors2 = np.sqrt(np.diag(cov2))

interc2 = unp.uarray(params2[0], errors2[0])
slope2 = unp.uarray(params2[1], errors2[1])

W2 = slope2 * sc.k / sc.elementary_charge

Tm = T[np.argmax(Ik)]
t0 = tau0(Tm, W2 * sc.e,  b)

# plotten
plt.figure(1)  # T, I 1
plt.plot(T, I, 'rx', label=r'$\mathrm{unkorrigierte} \ \mathrm{Werte}$')
plt.plot(T, Ik, 'x',  label=r'$\mathrm{korrigierte} \ \mathrm{Werte}$')
plt.plot(Tc, Ic, 'gx', label=r'$\mathrm{Werte} \ \mathrm{fuer} \ \mathrm{Anfangskurve}$')
plt.xlabel(r'$T / \mathrm{K}$')
plt.ylabel(r'$I / \mathrm{A}$')


plt.figure(2)  # Anfangsfit
x_plot = np.linspace(np.min(Tc), np.max(Tc), 1000)
plt.plot(1/Tc, np.log(jc), 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln(j \cdot \mathrm{m^2} \ / \ \mathrm{A}) $')
# plt.plot(1/T, np.log(I) - ground, 'gx')
plt.plot(1/x_plot, lin(1/x_plot, *params), 'b', label=r'$\mathrm{linearer} \ \mathrm{Fit}$')
plt.legend(loc='best')
plt.savefig("plots/lin1.pdf")

plt.figure(1)
ux_plot = np.linspace(np.min(Tu), np.max(Tu), 1000)
plt.plot(ux_plot, lin(ux_plot, *uparams), 'm', label=r'$\mathrm{Korrekturfit}$')
plt.legend(loc='best')
plt.savefig('plots/ti1.pdf')

plt.figure(3)
x_plot = np.linspace(np.min(T[mingrenz:maxgrenz]), np.max(T[mingrenz:maxgrenz]), 1000)
plt.plot(1/T, np.log(l), 'kx', label='umgerechnete Messwerte')
plt.plot(1/x_plot, lin(1/x_plot, *params2), 'b', label='linearer Fit')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln\frac{\mathrm{integral}}{I(T)}$')
plt.legend(loc='best')
plt.savefig('plots/lnlin1.pdf')

# ausgabe
# W mit methode 1
print("\nslope: ", slope)
print("\nintercept: ", interc)
print("\nW: ", W)

# W mit methode 2
print("\n\nslope2: ", slope2)
print("\nintercept2: ", interc2)
print("\nW2: ", W2)

print("\n\n\nintegral: ", integ)

print("\n\ntau0: ", t0)
print("Tmax: ", Tm)
print("b: ", b)
# daten abspeichern
np.savetxt('werte11.txt', np.column_stack([Ik[mingrenz:maxgrenz], T[mingrenz:maxgrenz]]), header="I korrigiert, T")
np.savetxt('werte12.txt', np.column_stack([Ic, Tc]), header="I anfang, T anfang")
# np.savetxt('werte13.txt', np.column_stack([]), header="")
# np.savetxt('werte14.txt', np.column_stack([,]), header="")
