import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties import unumpy as unp
import scipy.constants as sc
from scipy.optimize import curve_fit

# lade messwerte
spec = np.genfromtxt("lebenszeit.txt", unpack=True)


cal_spec = np.genfromtxt("kalibrierung.txt", unpack=True)
cal_peaks = cal_spec[cal_spec > 0]

cal_time = np.array(list(range(9)))
cal_time = (cal_time + 1.0) * 10**(-6)
print(cal_time)

tmess = 169098
# funktionen
def findpeak(spectrum):
    pre = 0
    peaks = []
    peaks_pos = []
    for i, val in enumerate(spectrum):
        if ((pre * val) != 0):
            peaks_pos.append((i * val / (val + pre) + (i - 1) *
                              pre / (val + pre)))
            peaks.append(val + pre)
            pre = val
        elif(val != 0 and spectrum[i + 1] == 0):
            peaks.append(val)
            peaks_pos.append(i)
        pre = val
    return(np.array(peaks_pos), np.array(peaks))


def lin(x, slope, intercept):
    return(x * slope + intercept)


def zerfall(t, n_0, lamb, U):
    return(n_0 * np.e**(- lamb * t) + U)


# rechnungen
# filtere peaks aus kalibrierungs-spectrum
spec_pos, cal_counts = findpeak(cal_spec)

# kalibrierung der zeitskala
params_cal, covariance_cal = curve_fit(lin, spec_pos, cal_time)
uparams_cal = unp.uarray(params_cal, np.sqrt(np.diag(covariance_cal)))

slope_cal = uparams_cal[0]
interc_cal = uparams_cal[1]

# wende kalibrierung an
real_time = lin(np.array(list(range(len(spec)))), slope_cal, interc_cal)

# fit der messwerte
params, covariance = curve_fit(zerfall, unp.nominal_values(real_time), spec, maxfev=2000)
uparams = unp.uarray(params, np.sqrt(np.diag(covariance)))

# plots
print("real t: ", real_time, "\ncal_peaks: ", cal_peaks)
# kalibrierung fit
plt.plot(spec_pos, cal_time, 'rx', label=r'$\mathrm{Messdaten für Kalibrierung}$')
plt.plot(spec_pos, lin(spec_pos, *params_cal), 'b', label=r'$\mathrm{linearer Fitfunktion}$')
plt.show()
plt.bar(unp.nominal_values(real_time), cal_spec, width=unp.nominal_values(slope_cal))
plt.show()


# fit an spectrum
plt.bar(unp.nominal_values(real_time), spec, width=unp.nominal_values(slope_cal))
# plt.plot(unp.nominal_values(real_time), spec)
plt.plot(unp.nominal_values(real_time), zerfall(unp.nominal_values(real_time), *params), 'r')
plt.show()

# ausgeben
print("spec: ", len(spec),  len(real_time))
print("fit: ", params)
# print("\nspec_pos: ", spec_pos, "\nspec_peaks: ", spec_pos)
# print("\n counts: ", cal_peaks)
