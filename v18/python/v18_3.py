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

X=np.array([121.78,244.70,295.94, 344.30, 411.12, 443.96, 678.00, 688.67, 778.90, 867.37, 964.08,1005.30,1085.90,1112.10,1299.10,1408.00,1457.60])

def f(x,a,b):
    return a*x+b

params,covariance=curve_fit(f,X,Y)

#Effizienzparameter
effizienz=unp.uarray(np.array([84.18399795734446,-1.028776078711911]), np.array([36.654658037116405, 0.07332701822586613]))
#36.654658037116405,0.07332701822586613
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

def intensitaet(peak_inhalt,peak_lage):
    return peak_inhalt/3238/(effizienz[0]*np.power(peak_lage,effizienz[1]))

# rechnungen
# b) analyse des photopeaks
'''
p2daten = gauss(daten2, 2220, 18)
halbwert = peakanteil(p2daten[0], 0.5)
zehntelwert = peakanteil(p2daten[0], 0.1)
halbwertth = halbwerttheorie(Ebild, p2daten[1])
zehntelwertth = zehntelwerttheorie(halbwertth)
p2inhalt = peakinhalt(p2daten)
comptonkante = comptkante(p2daten[1])
ruekstreupos = rueckstreu(p2daten[1])

minfit = 500
maxfit = 1570
comptparams, comptcovariance = curve_fit(kontinuum, np.linspace(minfit, maxfit, maxfit - minfit), daten2[minfit:maxfit])
ucomptparams = unp.uarray(comptparams, np.sqrt(np.diag(comptcovariance)))
# b) analyse des compton-kontinuums
# rstreupeak =
'''
# d) 
#Lese ungefaehre Peaks ab
d_guess = np.array([185, 277.67, 544, 753, 931.34, 1019.90, 1197.55, 1290.65])

#Berechne genaue Daten fuer jeden Peak

d_peak_breite = unp.uarray(np.zeros(len(d_guess)),np.zeros(len(d_guess))) #p1
d_peak_lage = unp.uarray(np.zeros(len(d_guess)),np.zeros(len(d_guess)))   #p2
d_peak_hoehe = unp.uarray(np.zeros(len(d_guess)),np.zeros(len(d_guess)))  #p3

#Fitbreite der for schleife unten
s=np.array([10,20, 5, 5,20,20,20,20])
for i in range(len(d_guess)):
    d_peak_breite[i] = gauss(daten3, d_guess[i],s[i])[0]
    d_peak_lage[i] = gauss(daten3, d_guess[i],s[i])[1]
    d_peak_hoehe[i] = gauss(daten3, d_guess[i],s[i])[3]
print('(d) Lage',d_peak_lage)
print('(d) Peakbreite :',d_peak_breite)
print('(d) Peakhoehe :',d_peak_hoehe)
#Berechne die Peakinhalte 
d_inhalte = d_peak_hoehe * unp.sqrt(np.pi / d_peak_breite)
print('(d) Peakinhalte :',d_inhalte)



#Berechne die Lage der Peaks in eV
d_lage_ev=umrechnen(d_peak_lage)
#print('(d) Peaklage Kanaele',d_peak_lage)
#print('(d) Peaklage in eV:',d_lage_ev)

#Intensitaet berechnen
d_peak_int=intensitaet(d_inhalte,d_lage_ev)
print('Intensitaeten :',d_peak_int)

#Aktivitaetsbestimmung
wkeit=np.array([2.2,34.06,0.65,0.45,7.2,18.3,62.1,8.9])/100

omega=0.2

Aktivitaet=4*np.pi/omega*d_peak_int/wkeit

print('Aktivitaet :',Aktivitaet)

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
#plt.figure(3)
#plt.bar(list(range(2300)), daten3[:2300], color='r')
#plt.plot(list(range(1600)), kontinuum(list(range(1600)), comptparams))
#
# ergebnisse
'''
print("aufgabe b): \n\n")
print("photpeak2: \nposition: ", p2daten[1])
print("\nhalbwert: ", halbwert)
print("\n10tel-wert: ", zehntelwert)
print("\nhablwert theorie: ", halbwertth)
print("\nzehntelwert theorie: ", zehntelwertth)
print("\npeakinhalt: ", p2inhalt)
print("\nzehntel/halbwert: ", zehntelwert/halbwert)
'''

np.savetxt('aktivitaet.txt',np.array([unp.nominal_values(d_peak_lage),unp.std_devs(d_peak_lage),unp.nominal_values(d_peak_int),unp.std_devs(d_peak_int),wkeit*100,unp.nominal_values(Aktivitaet),unp.std_devs(Aktivitaet)]).T,delimiter=' & ',newline=' ;newline; ',fmt=['%.2f',' %+.2f ',' %.2f ',' %+.2f ',' %.2f ',' %2d ',' %+2d'])
A=np.mean(unp.nominal_values(Aktivitaet[2:]))
A_std=np.std(unp.nominal_values(Aktivitaet[2:]))
print('Aktivitaet Mittelwert',unp.uarray(A,A_std))