import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from uncertainties import unumpy as unp
from uncertainties import ufloat

# daten einlesen
tEu = ufloat(4943, 5)
tdiff = 6084

daten1 = np.genfromtxt("v18_1.txt", unpack=True)
daten2 = np.genfromtxt("v18_2.txt", unpack=True)
daten3 = np.genfromtxt("v18_3.txt", unpack=True)
daten4 = np.genfromtxt("v18_4.txt", unpack=True)

# rechnungen


# plotten
# a)
plt.bar(list(range(len(daten1))), daten1, color='b')
plt.show()
# ergebnisse
