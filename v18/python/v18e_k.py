import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit
from scipy.optimize import fsolve as fs
from scipy.optimize import minimize as mini
from uncertainties import unumpy as unp
from uncertainties import ufloat
import scipy.integrate as integ
from scipy.optimize import minimize as minim

# einlesen der messwerte
daten1 = np.genfromtxt("v18_1.txt", unpack=True)
daten2 = np.genfromtxt("v18_2.txt", unpack=True)
daten3 = np.genfromtxt("v18_3.txt", unpack=True)
daten4 = np.genfromtxt("v18_4.txt", unpack=True)


Ebild = 2.9 # bildungsenergie
epsilon = 1 / (sc.m_e * sc.c**2)  # reziproke ruheenergie electron
krpeak2 = ufloat(638, 10)  # pos(kanal) rueckstreupeak
kckante2 = ufloat(1564, 10)  # pos(kanal) compt-kante

peak_pos = np.array([1846,
                     1764,
                     1154, 805, 700, 665, 75])

peak_pos *= 10**(3)
# Kalibrierung
Y = np.array([414.085320642,
             825.610503916,
             996.656026607,
             1158.50557442,
             1382.07913025,
             1491.85090388,
             2273.6046823,
             2310.7129048,
             2612.40170034,
             2908.13276148,
             3231.86620663,
             3368.34965244,
             3639.14375894,
             3726.05188115,
             4353.12130571,
             4716.3497603,
             4888.96958224])


X = np.array([121.78, 244.70, 295.94, 344.30, 411.12, 443.96, 678.00, 688.67,
              778.90, 867.37, 964.08, 1005.30, 1085.90, 1112.10, 1299.10,
              1408.00, 1457.60])
X *= 10**3


def f(x, a, b):
    return(x * a + b)


def antif(y, a, b):
    return((y - b) / a)


params, covariance = curve_fit(f, Y, X)  # kalibrierungs-fit

# funktionen
# gauss-funktion
def g(x, p1, p2, p3, p4):
    return(p4 * np.e**(-p1 * (x - p2)**2) + p3)


# gauss-fit an peak
def gauss1(data, pos, bereich):
    X2 = np.linspace(pos - bereich, pos + bereich, 2 * bereich + 1)
    #print("x2 : ", X2)
    kanalpos = np.round(antif(pos, *params))
    rkanalpos = kanalpos.astype(int)
    Y2 = data[rkanalpos - bereich:rkanalpos + bereich + 1]
    #print("y2 : ", Y2)
    params2, covariance2 = curve_fit(g, X2, Y2, p0=[1/900000, pos, 0, 1500], maxfev=5000)
    # Plot Vorbereiten
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 16
    plt.plot(X2, Y2, 'rx', label='Messwerte')
    x_plot = np.linspace(pos - bereich, pos + bereich, 500)
    plt.plot(x_plot, g(x_plot, *params2), label='Gaussfit')
    plt.xlabel(r'$\mathrm{Energie} \ / \ \mathrm{eV}$')
    plt.ylabel('Zaehlergebnis')
    plt.legend(loc='best')
    plt.show()
    uparams2 = unp.uarray(params2, np.sqrt(np.diag(covariance2)))
    # plt.savefig('plots/auswertung_b1.pdf'
    # halb = Y2[np.min(np.abs(Y2 - max_val / 2))]
    # zehnt = Y2[np.min(np.abs(Y2 - max_val / 10))]
    return(uparams2)

# rechnungen
gpeak_pos = np.zeros(len(peak_pos))
gpeak_pos_err = np.zeros(len(peak_pos))
peak_interv = np.zeros(len(peak_pos))
peak_interv = peak_interv.astype(int)
peak_interv += 10

def setinterv(index, wert):
    peak_interv[index] = wert
    return()

print("peak intervalle: ", peak_interv[2])
# setinterv()

for j, val in enumerate(peak_pos):
    h = gauss1(daten4, val, peak_interv[j])
    print("\n", h[1])
    gpeak_pos[j] = unp.nominal_values(h[1])
    gpeak_pos_err[j] = unp.std_devs(h[1])
    #gpeak_pos = ufloat(h[1].n, h[1].s)
rpeak_pos = unp.uarray(gpeak_pos, gpeak_pos_err)
print(ugpeak_pos)
print("\n\nreal postions: ", rpeak_pos)
