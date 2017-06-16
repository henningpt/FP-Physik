import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit
from scipy.optimize import fsolve as fs
from scipy.optimize import minimize as mini
from uncertainties import unumpy as unp
from uncertainties import ufloat
import scipy.integrate as integ
# daten einlesen
# tEu = ufloat(4943, 5)
# tdiff = 6084
Ebild = 2.9
epsilon = 1 / (sc.m_e * sc.c**2)
krpeak2 = ufloat(638, 10)

kckante2 = ufloat(1564, 10)

daten1 = np.genfromtxt("v18_1.txt", unpack=True)
daten2 = np.genfromtxt("v18_2.txt", unpack=True)
daten3 = np.genfromtxt("v18_3.txt", unpack=True)
daten4 = np.genfromtxt("v18_4.txt", unpack=True)

#Kalibrierung

Y=np.array([414.085320642,
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

X = np.array([121.78, 244.70, 295.94, 344.30, 411.12, 443.96, 678.00, 688.67, 778.90, 867.37, 964.08, 1005.30, 1085.90, 1112.10, 1299.10, 1408.00, 1457.60])


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
    plt.plot(X2, Y2, 'rx', label='Messdaten')
    plt.xlabel('Kanaele')
    plt.ylabel('Zaehlergebnis')
    plt.plot(np.linspace(np.min(X2), np.max(X2), 1000), g(np.linspace(np.min(X2), np.max(X2), 1000), *params2), label='Gaussfunktion')
    plt.legend(loc='best')
    plt.savefig('plots/gaussfit.pdf')

    return(uparams2)


def peakanteil(p1, pct):
    return(2 * unp.sqrt(unp.log(1/pct) / p1))


def peakinhalt(p):
    return(p[3] * unp.sqrt(params[0] *  np.pi / p[0]))


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
rpeak2 = umrechnen(krpeak2)
ckante2 = umrechnen(kckante2)
# b) analyse des photopeaks

p2daten = gauss(daten2, 2220, 18)
halbwert = peakanteil(p2daten[0] , 0.5) * params[0]
zehntelwert = peakanteil(p2daten[0] , 0.1) * params[0]
halbwertth = halbwerttheorie(Ebild, umrechnen(p2daten[1]))
zehntelwertth = zehntelwerttheorie(halbwertth)
p2inhalt = peakinhalt(p2daten)
comptonkante = comptkante(sc.e * 1000 * umrechnen(p2daten[1])) / sc.e
ruekstreupos = rueckstreu(sc.e * 1000 * umrechnen(p2daten[1])) / sc.e

minfit = 800
maxfit = 1570
comptparams, comptcovariance = curve_fit(kontinuum, np.linspace(minfit, maxfit, maxfit - minfit), daten2[minfit:maxfit])
ucomptparams = unp.uarray(comptparams, np.sqrt(np.diag(comptcovariance)))
def kontint(E):
    klammer = 2 + (E/(unp.nominal_values(p2daten[1]) - E))**2 * (1 / (unp.nominal_values(p2daten[1]) * epsilon)**2 + (unp.nominal_values(p2daten[1]) - E) / unp.nominal_values(p2daten[1]) - 2 * (unp.nominal_values(p2daten[1]) - E)/(unp.nominal_values(p2daten[1]) * unp.nominal_values(p2daten[1]) * epsilon))
    vorfakt = np.pi * sc.e**2 / (4 * np.pi * sc.epsilon_0 * sc.c**2 * sc.m_e)**2
    return(comptparams[0] * vorfakt * klammer)

comptinhalt = integ.quad(kontint, 55, 1570)[0] / params[0]


# b) analyse des compton-kontinuums
# rstreupeak =

# plotten
# spektren
'''
plt.figure(1)
plt.bar(list(range(len(daten1))), daten1, color='b')
plt.xlabel('Kanaele')
plt.ylabel('Zaehlergebnis')
plt.xlim((0,5000))
plt.ylim((0,400))
plt.savefig("plots/spec1.pdf")

plt.figure(2)
plt.bar(list(range(len(daten2))), daten2, color='b')
plt.xlabel('Kanaele')
plt.ylabel('Zaehlergebnis')
plt.xlim((0,2500))
plt.ylim((0,200))
plt.savefig("plots/spec2.pdf")

plt.figure(3)
plt.bar(list(range(len(daten3))), daten3, color='b')
plt.xlabel('Kanaele')
plt.ylabel('Zaehlergebnis')
plt.xlim((0,1500))
plt.ylim((0,500))
plt.savefig("plots/spec3.pdf")


plt.figure(4)
plt.bar(list(range(len(daten4))), daten4, color='b')
plt.xlabel('Kanaele')
plt.ylabel('Zaehlergebnis')
plt.xlim((0,5000))
plt.ylim((0,1000))
plt.show()
'''
plt.figure(2)
plt.bar(list(range(2300)), daten2[:2300], color='k', label='Messdaten')
plt.plot(list(range(1600)), kontinuum(list(range(1600)), comptparams), label='Fit des Kontinuums')
plt.xlabel('Kanaele')
plt.ylabel('Zaehlergebnis')
plt.xlim((0,1800))
plt.ylim((0,200))
plt.legend(loc='best')
plt.savefig('plots/kontinuumfit.pdf')

# ergebnisse
print("aufgabe b): \n\n")
print("photpeak2: \nposition: ", umrechnen(p2daten[1]))
print("\nhalbwert: ", halbwert)
print("\n10tel-wert: ", zehntelwert)
print("\nhablwert theorie: ", halbwertth)
print("\nzehntelwert theorie: ", zehntelwertth)
print("\npeakinhalt: ", p2inhalt)
print("\nzehntel/halbwert: ", zehntelwert/halbwert)
print("\nhalbwert rel abweichung: ", abs(halbwert - halbwertth) / halbwertth)
print("\nzehntelwert rel abweichung: ", abs(zehntelwert - zehntelwertth) / zehntelwertth)
print("\n\ncompton-kontinuum: \ncomptonkante: ", comptonkante)
print("comptonkante abgelesen: ", ckante2)
print("rel abweichung compton kante: ", abs(ckante2 *1000 - comptonkante)/( comptonkante))
print("\nrueckstreupeak abgelesen: ", rpeak2)
print("\nrueckstreupeak: ", ruekstreupos)
print("relative abweichung rpeak: ", abs(rpeak2 * 1000 - ruekstreupos)/(ruekstreupos))
print("\ncomptoninhalt: ", comptinhalt)
print("\ncompton paramter k: ", ucomptparams[0])
print("\n\nverhaeltnis photo zu continuum: ", p2inhalt / comptinhalt)
print('Fertig')
