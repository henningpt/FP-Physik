import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit
from scipy.optimize import fsolve as fs
from scipy.optimize import minimize as mini
from uncertainties import unumpy as unp
from uncertainties import ufloat
import scipy.integrate as integ

# Daten einlesen
daten1 = np.genfromtxt("v18_1.txt", unpack=True)
daten2 = np.genfromtxt("v18_2.txt", unpack=True)
daten3 = np.genfromtxt("v18_3.txt", unpack=True)
daten4 = np.genfromtxt("v18_4.txt", unpack=True)

peak_pos = np.array([4717, 4615, 4292, 4151, 3754,
                     3355, 3132, 2813, 2703, 2635, 2575, 2234,
                     2044, 994, 816, 629, 265, 1183, 488, 317, 337])


# Kallibrierung
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
             778.90, 867.37, 964.08, 1005.30,
             1085.90, 1112.10, 1299.10, 1408.00, 1457.60])


def f(x, a, b):
    return(x * a + b)

params, covariance = curve_fit(f, X, Y)


# funktionen
def g(x, p1, p2, p3, p4):
    return p4 * np.e**(-p1 * (x - p2)**2) + p3


def umrechnen(kanal):
    return kanal/params[0]-params[1]/params[0]


def gauss(spektrum, a, s):
    X2 = np.linspace(a - s, a + s, 2 * s + 1)
    # print('X2',X2)
    Y2 = spektrum[a - s:a + s + 1]
    # print('Y2',Y2)
    params2, covariance2 = curve_fit(g, X2, Y2,  p0=[1, a, 0, 1])
    uparams2 = unp.uarray(params2, np.sqrt(np.diag(covariance2)))
    # Plot Vorbereiten

    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 16
    plt.plot(X2, Y2, 'rx')
    plt.plot(np.linspace(np.min(X2), np.max(X2), 1000),
             g(np.linspace(np.min(X2), np.max(X2), 1000), *params2))
    plt.show()

    return(uparams2)


# rechnungen
gpeak_pos = np.zeros(len(peak_pos))
gpeak_pos_err = np.zeros(len(peak_pos))
peak_interv = np.zeros(len(peak_pos))
peak_interv = peak_interv.astype(int)
peak_interv += 15
def setinterv(index, wert):
    peak_interv[index] = wert



setinterv(15, 9)
setinterv(20, 12)
setinterv(19, 7)
setinterv(18, 7)
setinterv(16, 6)

for j, val in enumerate(peak_pos):
    h = gauss(daten4, val, peak_interv[j])
    print("\n", h[1])
    gpeak_pos[j] = unp.nominal_values(h[1])
    gpeak_pos_err[j] = unp.std_devs(h[1])
    #gpeak_pos = ufloat(h[1].n, h[1].s)
ugpeak_pos = unp.uarray(gpeak_pos, gpeak_pos_err)
print(ugpeak_pos)
rpeak_pos = umrechnen(ugpeak_pos)
print("\n\nreal postions: ", rpeak_pos)
