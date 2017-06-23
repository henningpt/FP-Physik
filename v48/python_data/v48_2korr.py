import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import scipy.constants as sc
import uncertainties.unumpy as unp
import scipy.integrate as integrate
from uncertainties import ufloat


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

b = 1.3 / 60
# funktionen


def kreisf(r):
    return(np.pi * r**2)


def nonegatives(arr):
    l = len(arr)
    for i in range(0, l):
        if (arr[i] < 0):
            arr[i] = 0.1 * 10**(-12)
    return(arr)


def lin(Tinv, intercept, slope):
    return(intercept + Tinv * slope)


def tau0(Tmax, W, b):
    return(sc.k * Tmax**2 * np.e**(-W / (sc.k * Tmax)) / (W * b))


def efct(T, a, b):
    return(a * np.e**(T * b))


def heizrate(T, zeit):
    summe = 0
    vorgaen = 0
    for i, val in enumerate(T):
        if (i > 0):
            summe += (abs(val - vorgaen) / zeit)
            vorgaen = val
    return(summe / len(T))


# rechnung
A = kreisf(1.5*10**(-3))


Tinv = 1/Tsort
# Ic = np.delete(Isort[:20], 5)
# Tc = np.delete(Tsort[:20], 5)

Torg = Tsort
Iorg = Isort

Tsort = np.delete(Tsort, [34, 35, 36])
Isort = np.delete(Isort, [34, 35, 36])

# untergrundfit + abziehen
Iu = np.append(Isort[9:10], Isort[68:76])
Tu = np.append(Tsort[9:10], Tsort[68:76])

params_u, covariance_u = curve_fit(lin, Tu, Iu)
eparams_u = unp.uarray(params_u, np.sqrt(np.diag(covariance_u)))
slope_u = eparams_u[1]
interc_u = eparams_u[0]

Ik = Isort - lin(Tsort, *params_u)
Ik2 = Ik[Isort >= Ik]
Tk2 = Tsort[Isort >= Ik]
Ir = Ik2[Ik2 >= 0]
Tr = Tk2[Ik2 >= 0]

# Anfangsfit
lin_fits = 11
lin_fite = 28
Ic = Ir[lin_fits:lin_fite]
jc = Ic/A  # stromdichte
Tc = Tr[lin_fits:lin_fite]
params_lin, cov_lin = curve_fit(lin, 1/Tc, np.log(Ic))
errors_lin = np.sqrt(np.diag(cov_lin))
interc_lin = ufloat(params_lin[0], errors_lin[0])
slope_lin = ufloat(params_lin[1], errors_lin[1])
W = -slope_lin * sc.k / sc.elementary_charge

# integration
intg_s = 8
intg_e = 64
integ = np.zeros(len(Ir[intg_s:intg_e]))
for i, val in enumerate(Ir[intg_s:intg_e]):
    integ_a = intg_s + i
    integ[i] = integrate.trapz(Ir[integ_a:intg_e], Tr[integ_a:intg_e]) / val

integ = np.delete(integ, len(integ) - 1)
ln_integ = np.log(integ)

params_int, covariance_int = curve_fit(lin, 1/Tr[intg_s:intg_e-1], ln_integ)
uparams_int = unp.uarray(params_int, np.sqrt(np.diag(covariance_int)))
slope_int = uparams_int[1]
interc_int = uparams_int[0]

W2 = slope_int * sc.k / sc.elementary_charge

# tau
bnew1 = heizrate(Tr[:15], 60)
bnew2 = heizrate(Tr[16:59], 30)
bnew3 = heizrate(Tr[60:], 60)
bnew = (bnew1 * len(Tr[:15]) + bnew2 * len(Tr[16:59]) + bnew3 * len(Tr[60:])) / len(Tr)
print("bnew: ", bnew)
Tm = Tr[np.argmax(Ir)]
t0 = tau0(Tm, W2 * sc.e,  bnew)


# plotten
# werte auswahl etc.
plt.figure(1)  # T, I 1
plt.plot(Torg, Iorg, 'yx', label=r'$\mathrm{herausgenommene} \ \mathrm{Messwerte}$')
plt.plot(Tsort, Isort, 'rx', label=r'$\mathrm{unkorr.} \ \mathrm{Werte}$')
plt.plot(Tsort, lin(Tsort, *params_u), label=r'$\mathrm{Untergrund-Fit}$')
plt.plot(Tsort, Ik, 'gx',  label=r'$\mathrm{korr.} \ \mathrm{Werte}$')
plt.plot(Tr, Ir, 'x', color=(1.0, 0.0, 1.0),
         label=r'$\mathrm{pos.korr.} \ \mathrm{Werte}$')


# plt.plot(Tc, Ic, 'gx', label=r'$\mathrm{Werte} \ \mathrm{fuer} \ \mathrm{Anfangskurve}$')
plt.xlabel(r'$T / \mathrm{K}$')
plt.ylabel(r'$I / \mathrm{A}$')
plt.legend(loc='best')
plt.show()

# anfangsplot
plt.figure(2)
x_plot = np.linspace(np.min(Tc), np.max(Tc), 1000)
plt.plot(1/Tc, np.log(Ic), 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
plt.plot(1/x_plot, lin(1/x_plot, *params_lin), 'b', label=r'$\mathrm{linearer} \ \mathrm{Fit}$')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln(j \cdot \mathrm{m^2} \ / \ \mathrm{A}) $')
plt.legend(loc='best')
plt.show()

# plot lin fit integration
plt.figure(3)
x_plot = np.linspace(Tr[intg_s], Tr[intg_e - 1], 1000)
plt.plot(1/Tr[intg_s:intg_e-1], np.log(integ), 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
# plt.plot(1/Tr[intg_s:intg_e-1], integ, 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
plt.plot(1/x_plot, lin(1/x_plot, *params_int), 'b', label=r'$\mathrm{linearer} \ \mathrm{Fit}$')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln\frac{\mathrm{integral}}{I(T)}$')
plt.legend(loc='best')
plt.show()


# ausgabe
#lin fit ergebnis W
print("Aktivierungsenergie W methode1: ", W)
print("\nAktivierungsenergie W methode2: ", W2)
print("tau0: ", t0)
# print("Ir: ", Ir)
# print("\nTr: ", Tr - 273.15)
