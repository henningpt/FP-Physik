import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import curve_fit
# from uncertainties import unumpy as unp
# from uncertainties import ufloat

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
    # Plot Vorbereiten
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 16
    plt.plot(X2, Y2, 'rx')
    plt.plot(X2, g(X2, *params2))
    plt.show()
    return(params2)

# rechnungen
# analyse des photopeaks
p2daten = (gauss(daten2, 2220, 18))

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
