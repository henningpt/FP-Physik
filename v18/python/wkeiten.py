import numpy as np


# werte
mu_ph = 39
mu_co = 0.8
D = 3.9 * 10**(-2)

# funktionen


def wkeit(mu, d):
    return(1 - np.e**(-mu * d))

# rechnungen
w_ph = wkeit(mu_ph, D)
w_co = wkeit(mu_co, D)


# ausgeben
print("Wkeit Ph: ", w_ph)
print("Wkeit Co: ", w_co)
print("verhaeltnis: ", w_ph / w_co)
