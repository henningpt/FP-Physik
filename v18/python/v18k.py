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


params, covariance = curve_fit(f, Y, X)  # kalibrierungs-fit


# funktionen
# gauss-funktion
def g(x, p1, p2, p3, p4):
    return(p4 * np.e**(-p1 * (x - p2)**2) + p3)


# gauss-fit an peak
def gauss1(data, pos, bereich):
    X2 = f(np.linspace(pos - bereich, pos + bereich, 2 * bereich + 1), *params)
    Y2 = data[pos - bereich:pos + bereich + 1]
    realpos = f(pos, *params)
    print("position: ", np.round(realpos))
    params2, covariance2 = curve_fit(g, X2, Y2, p0=[1/900000, realpos, 0, 1500])
    # Plot Vorbereiten
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 16
    plt.plot(X2, Y2, 'rx', label='Messwerte')
    x_plot = f(np.linspace(pos - bereich, pos + bereich, 500), *params)
    plt.plot(x_plot, g(x_plot, *params2), label='Gaussfit')
    plt.xlabel(r'$\mathrm{Energie} \ / \ \mathrm{eV}$')
    plt.ylabel('Zaehlergebnis')
    plt.legend()
    print("1")
    plt.show()
    # plt.savefig('plots/auswertung_b1.pdf')
    print("2")
    max_val = np.max(Y2)
    def gauss_halb(x):
        return(abs(params2[3] * np.e**(-params2[0] * (x - params2[1])**2) + params2[2] - (max_val / 2)))
    def gauss_zehnt(x):
        return(abs(params2[3] * np.e**(-params2[0] * (x - params2[1])**2) + params2[2] - (max_val / 10)))
    x_h_links = minim(gauss_halb, x0=660700)
    print("halb links: ", x_h_links.x)
    x_h_rechts = minim(gauss_halb, x0=662800)
    print("halb rechts: ", x_h_rechts.x)
    x_z_links = minim(gauss_zehnt, x0=659800)
    print("zehnt links: ", x_z_links.x)
    x_z_rechts = minim(gauss_zehnt, x0=663700)
    print("zehnt rechts: ", x_z_rechts.x)
    # halb = Y2[np.min(np.abs(Y2 - max_val / 2))]
    # zehnt = Y2[np.min(np.abs(Y2 - max_val / 10))]
    print("\nmax value : ", max_val)

    return(params2, np.sqrt(np.diag(covariance2)), abs(x_h_links.x - x_h_rechts.x), abs(x_z_links.x - x_z_rechts.x))


# funktion zur berechnung des peakinhalts
def peakinhalt(fit, pos, bereich):
    count_data = f(np.linspace(pos - bereich, pos + bereich, 2 * bereich + 1), *params)
    print(count_data)
    inhalt = np.sum(g(count_data, *fit))
    return(inhalt)

def peakinhalt_f(p):
    return(p[3] * unp.sqrt(params[0] * np.pi / p[0]))


# theoriewert für die halbwertsbreite
def halbwerttheorie(Eel, Eg):
    return(2.35 * unp.sqrt(0.1 * Eel * Eg))


# zusammenhang zwischen halbwertsbreite und zehntelwertsbreite
def zehntelwertth(halbwert):
    return(1.823 * halbwert)


# ermittelt die pos der compon-kante
def comptkante(Eg):
    return(Eg * 2 * Eg * epsilon / (1 + 2 * Eg * epsilon))


# ermittelt die pos des rueckstreupeaks
def rueckstreu(Eg):
    return(Eg / (1 + Eg * epsilon * 2))

# ermittelt inhalt des compton Kontinuums
def ck_inhalt(fit, ming, maxg):
    energien = f(np.linspace(ming, maxg, maxg - ming + 1), *params)
    inhalt = np.sum(kontinuum(energien, *fit))
    return(inhalt)

# funktion zur beschreibung des compton-kontinuums
def kontinuum(E, dsig):
    klammer = 2 + (E/(unp.nominal_values(p2daten[1]) - E))**2 * (1 / (unp.nominal_values(p2daten[1]) * epsilon)**2 + (unp.nominal_values(p2daten[1]) - E) / unp.nominal_values(p2daten[1]) - 2 * (unp.nominal_values(p2daten[1]) - E)/(unp.nominal_values(p2daten[1]) * unp.nominal_values(p2daten[1]) * epsilon))
    vorfakt = np.pi * sc.e**2 / (4 * np.pi * sc.epsilon_0 * sc.c**2 * sc.m_e)**2
    return(dsig * vorfakt * klammer)


# funktion zur bestimmung der n-tel-wertsbreite
def peakanteil(p1, pct):
    return(2 * unp.sqrt(unp.log(1/pct) / p1))

# funktion zur bestimmung der ntel wertsbreite (messen)

# rechnungen
rpeak2 = f(krpeak2, *params)
ckante2 = f(kckante2, *params)


# b) analyse des photopeaks
p2daten, p2err, xhalb, xzehnt = gauss1(daten2, 2220, 14)
up2daten = unp.uarray(p2daten, p2err)

halbwertth = halbwerttheorie(Ebild, up2daten)
zehntelwertth = zehntelwertth(halbwertth)

p2_inhalt = peakinhalt(up2daten, 2220, 20)  # inhalt photopeak


# analyse compton kontinuum
comptonkante = comptkante(sc.e * p2daten[1]) / sc.e
ruekstreupos = rueckstreu(sc.e * p2daten[1]) / sc.e

minfit = 800
maxfit = 1570
comptparams, comptcovariance = curve_fit(kontinuum, f(np.linspace(minfit, maxfit, maxfit - minfit), *params), daten2[minfit:maxfit])
ucomptparams = unp.uarray(comptparams, np.sqrt(np.diag(comptcovariance)))

ck2_inhalt = ck_inhalt(ucomptparams, minfit, maxfit)


# theorie halbwertsbreite etc.
halb_th = halbwerttheorie(Ebild, up2daten[1])


# plotten
'''
plt.figure(2)
plt.bar(left=f(np.linspace(0, 2000, 2001), *params), height=daten2[0:2001], color='b', label='Messdaten', width= p2daten[1] / len(daten2[:2001]))
plt.plot(f(np.linspace(0, 1570, 1570 + 1), *params), kontinuum(f(np.linspace(0, 1570, 1570 + 1), *params), comptparams), label='Fit des Kontinuums', color='y')
plt.plot(f(np.linspace(800, 1570, 1570 - 800 + 1), *params), kontinuum(f(np.linspace(800, 1570, 1570 - 800 + 1), *params), comptparams), label='Fit des Kontinuums im Bereich dafür verwendeten Werte', color='r')
plt.xlabel(r'$\mathrm{Energien} \ / \ \mathrm{eV}$')
plt.ylabel('Zaehlergebnis')
# plt.xlim((0,1800))
# plt.ylim((0,200))
leg = plt.legend(loc='best')
leg.get_frame().set_alpha(0.0)
plt.savefig('plots/kontinuumfit.pdf')
'''
# comptonk , rpos abgelesen
r_pos_ab = ufloat(188900, 3000)
c_k_ab = ufloat(465600, 3000)
# ausgeben
th_verh = ufloat(25, 31)
print("\np2 inhalt: ", p2_inhalt)
print("\nphoto fit parameter: ", up2daten)
print("\np2 ck_inhalt: ", ck2_inhalt)
print("\ncompton param: ", ucomptparams)
verhaeltnis = ck2_inhalt / p2_inhalt
print("\n peak / cont: ", verhaeltnis)
print("\n rel abweichung verhaeltnis: ", abs(verhaeltnis - th_verh) / th_verh )
print("xhalb: ", xhalb)
print("xzehnt: ", xzehnt)
print("habl th: ", halb_th)
print("rel abweichung: ", abs(xhalb - halb_th) / halb_th)

print("\n\n rueckstreu: ", ruekstreupos)
print("Abweichung rpos: ", abs(ruekstreupos - r_pos_ab) / ruekstreupos)
print("\ncomptonkante: ", comptonkante)
print("ABweichung compton k: ", abs(comptonkante - c_k_ab) / comptonkante)
