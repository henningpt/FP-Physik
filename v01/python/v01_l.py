import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties import unumpy as unp
import scipy.constants as sc
# from scipy.stats import curve_Fit

# lade messwerte
spec = np.genfromtxt("lebenszeit.txt", unpack=True)


cal_spec = np.genfromtxt("kalibrierung.txt", unpack=True)
c_peaks = cal_spec[cal_spec > 0]

c_time = np.array([list(range(9))])
c_time = (c_time + 1.0) * 10**(-6)

# funktionen
def findpeak(spectrum):
    pre = 0
    for i, val in enumerate(spectrum):
        if (pre * val > 0):
            
