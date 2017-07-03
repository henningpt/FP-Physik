import numpy as np
import matplotlib.pyplot as plt


# werte einlesen
vzt, counts = np.genfromtxt("werte_vz.txt", unpack=True)

     

# werte plotten
plt.errorbar(vzt, counts, yerr =  np.sqrt(counts), fmt='bx')
plt.plot(vzt, counts, 'x', color='r')
plt.show()
