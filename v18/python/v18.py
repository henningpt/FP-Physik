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
Ebild = 2.9 * sc.e
epsilon = 1 / (sc.m_e * sc.c**2)
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
    uparams2 = unp.uarray(params2, np.sqrt(np.diag(covariance2)))
    # Plot Vorbereiten
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 16
    plt.plot(X2, Y2, 'rx')
    plt.plot(np.linspace(np.min(X2), np.max(X2), 1000), g(np.linspace(np.min(X2), np.max(X2), 1000), *params2))
    # plt.show()
    return(uparams2)


def peakanteil(p1, pct):
    return(2 * unp.sqrt(unp.log(1/pct)/p1))


def peakinhalt(p):
    return(p[3] * unp.sqrt(np.pi / p[0]))


def halbwerttheorie(Eel, Eg):
    return(2.35 * unp.sqrt(0.1 * Eel * Eg))


def zehntelwerttheorie(halbwert):
    return(1.823 * halbwert)


def comptkante(Eg):
    return(Eg * 2 * Eg * epsilon / (1 + 2 * Eg * epsilon))


def rueckstreu(Eg):
    return(Eg / (1 + Eg * epsilon * 2))


def kontinuum(E, dsig):
    klammer = 2 + (E/(unp.nominal_values(p2daten[1]) - E))**2 * (1 / (unp.nominal_values(p2daten[1]) * epsilon)**2 + (unp.nominal_values(p2daten[1]) - E) / unp.nominal_values(p2daten[1]) - 2 * (unp.nominal_values(p2daten[1]) - E)/(unp.nominal_values(p2daten[1]) * unp.nominal_values(p2daten[1]) * epsilon))
    vorfakt = np.pi * sc.e**2 / (4 * np.pi * sc.epsilon_0 * sc.c**2 * sc.m_e)**2
    return(dsig * vorfakt * klammer)


# rechnungen
# b) analyse des photopeaks
p2daten = gauss(daten2, 2220, 18)
halbwert = peakanteil(p2daten[0], 0.5)
zehntelwert = peakanteil(p2daten[0], 0.1)
halbwertth = halbwerttheorie(Ebild, p2daten[1])
zehntelwertth = zehntelwerttheorie(halbwertth)
p2inhalt = peakinhalt(p2daten)
comptonkante = comptkante(p2daten[1])
ruekstreupos = rueckstreu(p2daten[1])

minfit = 1000
maxfit = 2000
comptparams, comptcovariance = curve_fit(kontinuum, np.linspace(minfit, maxfit, maxfit - minfit), daten2[minfit:maxfit])

# b) analyse des compton-kontinuums
# rstreupeak =

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
# plt.figure(2)
# plt.bar(list(range(len(daten2))), daten2, color='r')
# plt.show()

# ergebnisse
print("aufgabe b): \n\n")
print("photpeak2: \nposition: ", p2daten[1])
print("\nhalbwert: ", halbwert)
print("\n10tel-wert: ", zehntelwert)
print("\nhablwert theorie: ", halbwertth)
print("\nzehntelwert theorie: ", zehntelwertth)
print("\npeakinhalt: ", p2inhalt)
print("\nzehntel/halbwert: ", zehntelwert/halbwert)
