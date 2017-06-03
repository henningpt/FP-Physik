import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit
from scipy.optimize import fsolve as fs
from scipy.optimize import minimize as mini
from uncertainties import unumpy as unp
from uncertainties import ufloat

# daten einlesen
# tEu = ufloat(4943, 5)
# tdiff = 6084

daten1 = np.genfromtxt("v18_1.txt", unpack=True)
daten2 = np.genfromtxt("v18_2.txt", unpack=True)
daten3 = np.genfromtxt("v18_3.txt", unpack=True)
daten4 = np.genfromtxt("v18_4.txt", unpack=True)


# funktionen
def g(x, p1, p2, p3, p4):
    return p4 * np.e**(-p1 * (x - p2)**2) + p3


def gauss(spektrum, a, s):
    X2 = np.linspace(a - s, a + s, 2 * s + 1)
    # print('X2',X2)
    Y2 = spektrum[a - s:a + s + 1]
    # print('Y2',Y2)
    params2, covariance2 = curve_fit(g, X2, Y2,  p0=[1, a, 0, 1])
    uparams2 = unp.uarray(params2,np.sqrt(np.diag(covariance2)))
    # Plot Vorbereiten
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 16
    plt.plot(X2, Y2, 'rx')
    plt.plot(np.linspace(np.min(X2), np.max(X2), 1000), g(np.linspace(np.min(X2), np.max(X2), 1000), *params2))
    plt.show()
    return(params2, uparams2)


def peakanteil(fun, params, pct):
    def usefun(x):
        return(params[3] * np.e**(-params[0] * (x - params[1])**2) + params[2])

    def minfun(x):
        return(pct * usefun(params[1]) - usefun(x))
    anteilv = fs(minfun, x0=params[1]*0.999)
    return(2*abs(params[1] - anteilv))


# rechnungen
# b) analyse des photopeaks
p2daten, p2test = gauss(daten2, 2220, 18)
halbwert = peakanteil(g, p2daten, 0.5)
zehntelwert = peakanteil(g, p2daten, 0.1)
# peakinhalt =
# plotten
# spektren
'''
plt.figure(1)
plt.bar(list(range(len(daten1))), daten1, color='r')
plt.savefig("plots/spec1.pdf")

plt.figure(2)
plt.bar(list(range(len(daten2))), daten2, color='r')
plt.savefig("plots/spec2.pdf")

plt.figure(3)
plt.bar(list(range(len(daten3))), daten3, color='r')
plt.savefig("plots/spec3.pdf")

plt.figure(4)
plt.bar(list(range(len(daten4))), daten4, color='r')
plt.savefig("plots/spec4.pdf")
'''

# ergebnisse
print("aufgabe b): \n\n")
print("photpeak2: \n position: ", p2daten[1])
print("halbwert: ", halbwert)
print("10tel-wert: ", zehntelwert)
print("test: ", p2test)
