import numpy as np
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants 
# daten einlesen
U_unc = 2  # fehler auf spannung +-2

um_0, U_0 = np.genfromtxt("m_0.txt", unpack=True)
um_1, U_1 = np.genfromtxt("m_1.txt", unpack=True)
um_2, U_2 = np.genfromtxt("m_2.txt", unpack=True)
um_3, U_3 = np.genfromtxt("m_3.txt", unpack=True)
um_4, U_4 = np.genfromtxt("m_4.txt", unpack=True)
um_5, U_5 = np.genfromtxt("m_5.txt", unpack=True)
um_6, U_6 = np.genfromtxt("m_6.txt", unpack=True)
um_7, U_7 = np.genfromtxt("m_7.txt", unpack=True)


# werte umrechnen
def umrechnen(umdreh):
    return(umdreh * 1.8 * 10**(-3))

s_0 = umrechnen(um_0)
s_1 = umrechnen(um_1)
s_2 = umrechnen(um_2)
s_3 = umrechnen(um_3)
s_4 = umrechnen(um_4)
s_5 = umrechnen(um_5)
s_6 = umrechnen(um_6)
s_7 = umrechnen(um_7)


# funktionen
def plot_SU(s, U, nummer):
    nstr = str(nummer)
    plt.figure(nummer)
    plt.plot(s, U, 'bx')
    plt.xlabel(r'$ s \ / \ \mathrm{m} $')
    plt.ylabel(r'$ U \ / \ \mathrm{V} $')
    plt.show()
   # plt.savefig("plots/plot_" + nstr)
    
def linfit(x,a,b):
    return a*x + b

# plotten
'''
plot_SU(s_0, U_0, 0)
plot_SU(s_1, U_1, 1)
plot_SU(s_2, U_2, 2)
plot_SU(s_3, U_3, 3)
plot_SU(s_4, U_4, 4)
plot_SU(s_5, U_5, 5)
plot_SU(s_6, U_6, 6)
plot_SU(s_7, U_7, 7)
'''
#Abgelesene Werte

null=0.0116306
max2=np.array([0.0123886,0.0125965,0.0127545,0.0129138,0.0131821,0.0132806])
max1=np.array([0.0107086,0.0105543,0.0103943,0.0102490,0.0100221,0.0099415])
wert1=null-max1
wert2=max2-null
Werte=(max2-max1)/2
#Werte=np.array([0.012383-0.010723,0.012540-0.010545,0.012737-0.0103925,0.012948-0.010251, 0.013183 - 0.009968, 0.013327 - 0.009948])/2

#Fit
a=0.0025
L=0.07
l=0.455
B=np.array([0.23,0.30,0.38,0.45, 0.59, 0.65])
Gradient=0.968*B/a
params, covariance = curve_fit( linfit, Gradient, Werte )  

uparams=unp.uarray(params,np.sqrt(np.diag(covariance)))

print('Steigung: ', uparams[0])

X=np.linspace(Gradient[0],Gradient[len(Gradient)-1],100)

plt.figure(8)
plt.plot(Gradient,Werte,'bx',label='Messwerte')
plt.plot(X, linfit(X,*params), 'r-',label='lineare Regression')
plt.xlabel(r'$ \frac{\partial B}{\partial z} \ /  \frac{\mathrm{T}}{\mathrm{m}} $')
plt.ylabel(r'$ s \ / \ \mathrm{m} $')
plt.legend()
#plt.show()
plt.savefig("plots/plot_8")

T=190+273.15
k=scipy.constants.Boltzmann
g=2
m=0.5

#Bohrsches Magneton
mu_berechnet=uparams[0]*6*k*T/(l*L*(1-L/(2*l)))/(g*m)
print('Bohrsches Magneton berechnet: ', mu_berechnet)

mu=scipy.constants.physical_constants['Bohr magneton']
print('Bohrsches Magneton Literatur: ', mu)

tab=np.array([B,max1*1000,max2*1000,Werte*1000]).T
np.savetxt('tab.txt', np.column_stack(tab).T,delimiter='&',fmt='%2f',newline=' xXx')

print('Abweichung: ',(mu[0]-mu_berechnet)/mu[0])

print('done.')

