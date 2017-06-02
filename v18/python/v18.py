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


# rechnungen


# plotten
# spektren
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


# ergebnisse
