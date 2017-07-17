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
# print(cal_time)

cstart = 3920300
tmess = 169098
tsuch = 20 * 10**(-6)

lifetime_lit = 2.19704 * 10**(-6)


# funktionen
def findpeak(spectrum):
    pre = 0
    peaks = []
    peaks_pos = []
    for i, val in enumerate(spectrum):
        if ((pre * val) != 0):
            peaks_pos.append(i * val / (val + pre) + (i - 1) *
                              pre / (val + pre))
            print("\nstandard error: ", np.std([i, i - 1]) * np.sqrt((val/(val+pre))**2 + (pre/(val+pre))**2) / np.sqrt(2))
            peaks.append(ufloat(val + pre, np.std([val, pre]) / np.sqrt(2)))
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


def zerfall2(t, n_0, lamb):
    return(n_0 * np.e**(- lamb * t))


def untergrund(c_start, t_mess, t_such):
    rate = c_start / t_mess
    erwartung = rate * t_such
    wkeit = erwartung * np.e**(-erwartung)
    c_untergrund = c_start * wkeit
    return(c_untergrund / 512)


def cut_zeros(arr):
    vorg = -1
    for i, val in enumerate(arr):
        if (vorg == 0 and val == 0):
            return(arr[:i])
    return(arr)


# rechnungen
# filtere peaks aus kalibrierungs-spectrum
spec_pos, cal_counts = findpeak(cal_spec)

# kalibrierung der zeitskala
params_cal, covariance_cal = curve_fit(lin, spec_pos, cal_time)
uparams_cal = unp.uarray(params_cal, np.sqrt(np.diag(covariance_cal)))

slope_cal = uparams_cal[0]
interc_cal = uparams_cal[1]

# untergrund berechnen
U_rate_c = untergrund(ufloat(cstart, np.sqrt(cstart)), tmess, tsuch)
U_rate = unp.nominal_values(U_rate_c)

# wende kalibrierung an
real_time = lin(np.array(list(range(len(spec)))), slope_cal, interc_cal)
spec_k = spec[spec > 0]
real_time_k = real_time[spec > 0]


# fit der messwerte
params, covariance = curve_fit(zerfall, unp.nominal_values(real_time_k),
                               spec_k, maxfev=800)
uparams = unp.uarray(params, np.sqrt(np.diag(covariance)))


params2, covariance2 = curve_fit(zerfall2,
                                 unp.nominal_values(real_time_k), spec_k - U_rate)

uparams2 = unp.uarray(params2, np.sqrt(np.diag(covariance2)))


# plots
# kalibrierung fit
plt.figure(1)
plt.plot(spec_pos, cal_time, 'rx',
         label=r'$\mathrm{Messdaten \ für \ Kalibrierung}$')
plt.plot(spec_pos, lin(spec_pos, *params_cal), 'b',
         label=r'$\mathrm{lineare \ Fitfunktion}$')
plt.ylabel(r'$t \ / \ \mathrm{s}$')
plt.xlabel(r'$\mathrm{Kanal}$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc='best')
plt.savefig("plots/kalibrierung.pdf")
plt.bar(unp.nominal_values(real_time), cal_spec,
        width=unp.nominal_values(slope_cal))



# fit an spectrum (untergrund berechnet)
# plt.bar(unp.nominal_values(real_time), spec - U_rate,
#        width=unp.nominal_values(slope_cal), color='b', edgecolor='none',
#        log=0, alpha=0.33, label=r'$\mathrm{Messdaten}$')
# plt.plot(unp.nominal_values(real_time), spec, 'b')
plt.figure(2)
plt.errorbar(unp.nominal_values(real_time), spec - U_rate, yerr=np.sqrt(spec - U_rate), fmt='x',
             color='r', alpha=0.4, label=r'$\mathrm{Messdaten - Untergrund}$', zorder=1)
plt.plot(unp.nominal_values(real_time_k),
         zerfall2(unp.nominal_values(real_time_k), *params2), 'g',
         label=r'$\mathrm{Exponentieller-Fit}$', alpha=1.0, zorder=2)
plt.xlabel(r'$t \ / \ \mathrm{s} $')
plt.ylabel(r'$\mathrm{counts}$')
plt.legend(loc='best')
# plt.yscale('log')
plt.savefig("plots/hist_u.pdf")

plt.figure(4)
plt.errorbar(unp.nominal_values(real_time), spec - U_rate, yerr=np.sqrt(spec - U_rate), fmt='x',
             color='r', alpha=0.4, label=r'$\mathrm{Messdaten - Untergrund}$', zorder=1)
plt.plot(unp.nominal_values(real_time_k),
         zerfall2(unp.nominal_values(real_time_k), *params2), 'g',
         label=r'$\mathrm{Fit \ an \ Messdaten}$', alpha=1.0, zorder=2)
plt.xlabel(r'$t \ / \ \mathrm{s} $')
plt.ylabel(r'$\ln(\mathrm{counts})$')
plt.yscale('log')
plt.legend(loc='best')
# plt.yscale('log')
plt.savefig("plots/hist_u_log.pdf")

# fit an spectrum (untergrund als fit parameter)
# plt.bar(unp.nominal_values(real_time), spec,
#        width=unp.nominal_values(slope_cal), color='b', edgecolor='none',
#        log=0, alpha=0.33, label=r'$\mathrm{Messdaten}$')
# plt.plot(unp.nominal_values(real_time), spec, 'b')
plt.figure(3)
plt.errorbar(unp.nominal_values(real_time), spec, yerr=np.sqrt(spec), fmt='x',
             color='r', alpha=0.4, label=r'$\mathrm{Messdaten}$', zorder=1)
plt.plot(unp.nominal_values(real_time_k),
         zerfall(unp.nominal_values(real_time_k), *params), 'g',
         label=r'$\mathrm{Exponentieller-Fit}$', alpha=1.0, zorder=2)
plt.xlabel(r'$t \ / \ \mathrm{s} $')
plt.ylabel(r'$\mathrm{counts}$')
plt.legend(loc='best')
# plt.yscale('log')
plt.savefig("plots/hist.pdf")


plt.figure(5)
plt.errorbar(unp.nominal_values(real_time), spec, yerr=np.sqrt(spec), fmt='x',
             color='r', alpha=0.4, label=r'$\mathrm{Messdaten}$', zorder=1)
plt.plot(unp.nominal_values(real_time_k),
         zerfall(unp.nominal_values(real_time_k), *params), 'g',
         label=r'$\mathrm{Fit \ an \ Messdaten}$', alpha=1.0, zorder=2)
plt.xlabel(r'$t \ / \ \mathrm{s} $')
plt.ylabel(r'$\ln(\mathrm{counts})$')
plt.legend(loc='best')
plt.yscale('log')
plt.savefig("plots/hist_log.pdf")

# ausgeben
print("calibration params: ", uparams_cal)
print("\nspec length: ", len(spec),  len(real_time))
print("\nfit: ", uparams)
print("fit2: ", uparams2)

# print("\nspec_pos: ", spec_pos, "\nspec_peaks: ", spec_pos)
# print("\n counts: ", cal_peaks)
print("\nUntergrund berechnet: ", U_rate_c)

print("\n\nlebensdauer 1: ", 1 / uparams[1])
print("rel. abweichung 1: ", abs(1 / uparams[1] - lifetime_lit) / lifetime_lit)
print("\nlebensdauer 2: ", 1 / uparams2[1])
print("rel. abweichung 2: ", abs(1 / uparams2[1] - lifetime_lit) / lifetime_lit)
print("\n\n spec_pos: ", spec_pos)


print("\n\nvalid kanaele: ", len(real_time_k))
print("gesamt counts: ", np.sum(unp.uarray(spec, np.sqrt(spec))))
