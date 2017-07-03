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

Tstart = 210.0
W_lit = 0.66
t_lit = 4 * 10**(-14)

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


def tau0(T_0, Tmax, W, b):
    return(sc.k * Tmax**2 * (np.e**(-W / (sc.k * Tmax)) -
           np.e**(-W / (sc.k * T_0))) / (W * b))


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

# untergrundfit + abziehen
Iu = np.append(I[10:12], I[28:30])
Tu = np.append(T[10:12], T[28:30])

params_u, covariance_u = curve_fit(lin, Tu, Iu)
eparams_u = unp.uarray(params_u, np.sqrt(np.diag(covariance_u)))
slope_u = eparams_u[1]
interc_u = eparams_u[0]

Ik = I - lin(T, *params_u)
Ik2 = Ik[I >= Ik]
Tk2 = T[I >= Ik]
Ir = Ik2[Ik2 >= 0]
Tr = Tk2[Ik2 >= 0]

# Anfangsfit
lin_fits = 4
lin_fite = 9
Ic = Ir[lin_fits:lin_fite]
jc = Ir/A  # stromdichte
Tc = Tr[lin_fits:lin_fite]
params_lin, cov_lin = curve_fit(lin, 1/Tc, np.log(Ic))
errors_lin = np.sqrt(np.diag(cov_lin))
interc_lin = ufloat(params_lin[0], errors_lin[0])
slope_lin = ufloat(params_lin[1], errors_lin[1])
W = -slope_lin * sc.k / sc.elementary_charge

deltaW1 = abs(W - W_lit) / W_lit

# integration
intg_s = 4
intg_e = len(Ir)
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

deltaW2 = abs(W2 - W_lit) / W_lit

# tau
bnew = heizrate(Tr, 60)
Tm = Tr[np.argmax(Ir)]
t0 = tau0(Tstart, Tm, W2 * sc.e,  bnew)

delta_t = abs(t0 - t_lit) / t_lit

# plotten
# werte auswahl etc.
plt.figure(1)  # T, I 1
plt.plot(T, I, 'rx', label=r'$\mathrm{unkorr.} \ \mathrm{Werte}$')
plt.plot(T, lin(T, *params_u), label=r'$\mathrm{Untergrund-Fit}$')
plt.plot(T, Ik, 'x', color=(1.0, 0.0, 1.0),  label=r'$\mathrm{korr.} \ \mathrm{Werte}$')
plt.plot(Tr, Ir, 'gx',
         label=r'$\mathrm{pos.korr.} \ \mathrm{Werte}$')

# plt.plot(Tc, Ic, 'gx', label=r'$\mathrm{Werte} \ \mathrm{fuer} \ \mathrm{Anfangskurve}$')
plt.xlabel(r'$T / \mathrm{K}$')
plt.ylabel(r'$I / \mathrm{A}$')
plt.legend(loc='best')
plt.savefig("plots/ti1.pdf")


# anfangsplot
plt.figure(2)
x_plot = np.linspace(np.min(Tc), np.max(Tc), 1000)
plt.plot(1/Tc, np.log(Ic), 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
plt.plot(1/x_plot, lin(1/x_plot, *params_lin), 'b', label=r'$\mathrm{linearer} \ \mathrm{Fit}$')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln(Ix  / \ \mathrm{A}) $')
plt.legend(loc='best')
plt.savefig("plots/lin1.pdf")

# plot lin fit integration
plt.figure(3)
x_plot = np.linspace(Tr[intg_s], Tr[intg_e - 1], 1000)
plt.plot(1/Tr[intg_s:intg_e-1], np.log(integ), 'rx', label=r'$\mathrm{umgerechnete} \ \mathrm{Messwerte}$')
plt.plot(1/x_plot, lin(1/x_plot, *params_int), 'b', label=r'$\mathrm{linearer} \ \mathrm{Fit}$')
plt.xlabel(r'$\frac{1}{T} \ / \ \frac{1}{\mathrm{K}}$')
plt.ylabel(r'$\ln\frac{\mathrm{integral}}{I(T)}$')
plt.legend(loc='best')
plt.savefig("plots/lnlin1.pdf")


# ausgabe
# lin fit ergebnis W
print("\n\nAktivierungsenergie W: ", W)
print("\nAktivierungsenergie W2: ", W2)
print("\n\ntau0: ", t0)
print("\nb: ", bnew)


print("\n\n Abweichung W1 : ", deltaW1)
print("\n Abweichung W2 : ", deltaW2)
print("\n Abweichung t0 : ", delta_t)

# paramter der Fits
print("\n\nslope untergrund : ", slope_u)
print("\ninterc untergrund : ", interc_u)

print("\n\nslope anfang : ", slope_lin)
print("\ninterc anfang : ", interc_lin)

print("\n\nslope integral : ", slope_int)
print("\ninterc integral : ", interc_int)

np.savetxt('werte_untergrund1.txt',np.column_stack([Iu, Tu]), header="I untergrund, T untergrund")

# werte fuer anfangsfit und integration abspeichern
np.savetxt('werte_anfang1.txt', np.column_stack([Ic, Tc]),
           header="I anfang, T anfang")

np.savetxt('werte_integration1.txt', np.column_stack([Ir[intg_s:intg_e], Tr[intg_s:intg_e]]), header="I anfang, T anfang")
# print("\nTr: ", Tr)
