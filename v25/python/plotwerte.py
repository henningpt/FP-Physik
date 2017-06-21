import numpy as np
import matplotlib.pyplot as plt

# daten einlesen
U_unc = 2  # fehler auf spannung +-2

um_0, U_0 = np.genfromtxt("m_0.txt", unpack=True)
um_1, U_1 = np.genfromtxt("m_1.txt", unpack=True)
um_2, U_2 = np.genfromtxt("m_2.txt", unpack=True)
um_3, U_3 = np.genfromtxt("m_3.txt", unpack=True)
um_4, U_4 = np.genfromtxt("m_4.txt", unpack=True)
um_5, U_5 = np.genfromtxt("m_5.txt", unpack=True)
# um_6, U_6 = np.genfromtxt("m_6.txt", unpack=True)
# um_7, U_7 = np.genfromtxt("m_7.txt", unpack=True)


# werte umrechnen
def umrechnen(umdreh):
    return(umdreh * 1.8 * 10**(-3))

s_0 = umrechnen(um_0)
s_1 = umrechnen(um_1)
s_2 = umrechnen(um_2)
s_3 = umrechnen(um_3)
s_4 = umrechnen(um_4)
s_5 = umrechnen(um_5)
# s_6 = umrechnen(um_6)
# s_7 = umrechnen(um_7)


# funktionen
def plot_SU(s, U, nummer):
    nstr = str(nummer)
    plt.figure(nummer)
    plt.plot(s, U, 'bx')
    plt.xlabel(r'$ s \ / \ \mathrm{m} $')
    plt.ylabel(r'$ U \ / \ \mathrm{V} $')
    plt.savefig("plots/plot_" + nstr)
# plotten
plot_SU(s_0, U_0, 0)
plot_SU(s_1, U_1, 1)
plot_SU(s_2, U_2, 2)
plot_SU(s_3, U_3, 3)
plot_SU(s_4, U_4, 4)
plot_SU(s_5, U_5, 5)
