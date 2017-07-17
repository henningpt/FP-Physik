import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.optimize import minimize as mini
from scipy.stats import sem
# werte einlesen
messz = 10  # seconds

vzt_chaos, counts_chaos = np.genfromtxt("werte_vz.txt", unpack=True)
# sortiere arrays
order = np.argsort(vzt_chaos)
vzt = 10**(-9) * vzt_chaos[order]
counts = counts_chaos[order] / messz
print("test : ", vzt)

# intervalle
vzt_1 = vzt[0:9]
vzt_2 = vzt[8:29]
vzt_3 = vzt[29:len(counts) + 1]

counts_1 = counts[0:9]
counts_2 = counts[8:29]
counts_3 = counts[29:len(counts) + 1]

# functions
def lin(x, a, b):
    return(x * a + b)


def const(x, c):
    return(c * x**(0))


# fits an messwerte

fit1, err1 = curve_fit(lin, vzt_1, counts_1)
fit2, err2 = curve_fit(const, vzt_2, counts_2)
fit3, err3 = curve_fit(lin, vzt_3, counts_3)

# finde maxwert + halbwert
umaxwert = ufloat(fit2[0], err2[0])
maxwert = umaxwert.n

def leftwert(x):
    return(lin(x, *fit1) - maxwert / 2)


def rightwert(x):
    return(lin(x, *fit3) - maxwert / 2)


lhalb = mini(leftwert, x0=-30.0 * 10**(-9)).x
rhalb = mini(rightwert, x0=30.0 * 10**(-9)).x

deltat =0.5 * (abs(lhalb) + abs(rhalb))

# werte plotten
plt.errorbar(vzt, counts, yerr=np.sqrt(counts / messz), fmt='bx', label=r'$\mathrm{Messwerte}$')
plt.plot(vzt, counts, 'x', color='r')
plt.plot(vzt_1, lin(vzt_1, *fit1), color='r', label=r'$\mathrm{Lineare \ Regression \ der \ Messwerte}$')
plt.plot(vzt_2, const(vzt_2, *fit2), color='g', label=r'$\mathrm{Konstante \ Regression \ der \ Messwerte}$')
plt.plot(vzt_3, lin(vzt_3, *fit3), color='c', label=r'$\mathrm{Lineare \ Regression \ der \ Messwerte}$')
plt.legend(loc='best')
plt.xlabel(r'$\Delta t \ / \ \mathrm{s}$')
plt.ylabel(r'$\mathrm{counts} \ / \ \mathrm{s}$')

plt.savefig("plots/verzoegerungszeit.pdf")


# ausgabe
print("maxwert: ", umaxwert)
print("\nleft halb: ", lhalb)
print("right halb: ", rhalb)
print("mittel habl: ", deltat)
