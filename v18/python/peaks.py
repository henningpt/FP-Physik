
import numpy as np
import math
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
WERTE=np.genfromtxt('v18_1.txt' , unpack=True)

Y=np.array([414.085320642,
825.610503916,
996.656026607,
1158.50557442,
1382.07913025,
1491.85090388,
2273.6046823,
2310.7129048,
2612.40170034,
2908.13276148,
3231.86620663,
3368.34965244,
3639.14375894,
3726.05188115,
4353.12130571,
4716.3497603,
4888.96958224])

X=np.array([121.78,244.70,295.94, 344.30, 411.12, 443.96, 678.00, 688.67, 778.90, 867.37, 964.08,1005.30,1085.90,1112.10,1299.10,1408.00,1457.60])

A=np.array([121,244,295,344,411,443,678,688,778,867,964,1005,1085,1112,1299,1408,1457])
#print(A)
def f(x,a,b):
    return a*x+b

params,covariance=curve_fit(f,X,Y)

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

plt.plot(np.linspace(100,1500,100),f(np.linspace(100,1500,100),*params),label='Fitfunktion')
plt.plot(X,Y,'rx',label='Messergebnisse')
plt.ylabel('Kanalnummer')
plt.xlabel('Energie in eV')
plt.legend()
#plt.show()

#Kalibrierungsparameter
kal=unp.uarray(params,np.sqrt(np.diag(covariance)))
print('Kalibrierungswerte :',kal)
#Funktionen
#Gauss
def gauss(a,s):
 X2=np.linspace(a-s,a+s,2*s+1)
# print('X2',X2)
 Y2=WERTE[a-s:a+s+1]
 #print('Y2',Y2)
 def g(x,p1,p2,p3,p4):
    return p4*np.e**(-p1*(x-p2)**2) + p3

 params2,covariance2 = curve_fit(g,X2,Y2, p0=[1,a,0,1])
 

 #Plot Vorbereiten
 import matplotlib.pyplot as plt
 #plt.rcParams['figure.figsize'] = (10, 8)
 #plt.rcParams['font.size'] = 16
 
 #plt.plot(X2,Y2,'rx')
 #plt.plot(X2,g(X2,*params2))
 #plt.show() 
 return params2, covariance2

#Abgelesene Y Werte mit Gauss fitten und in Y eintragen
C=[]
for i in range(len(Y)):
    c,d=gauss(Y[i],20)
    C=C+[c]
    Y[i]=C[i][1] 
#Fuer bessere Werte nochmal
C=[]
D=[]
for i in range(len(Y)):
    c,d=gauss(Y[i],20)
    C=C+[c]
    D=D+[d]
    Y[i]=C[i][1] 
#Zeige werte an
#for i in range(len(C)):
 #   print(Y[i])
    
#print('Parameter',gauss(413,20))

#Setzte Y Fehler
Y_err=np.zeros(len(Y))
for i in range(len(Y)):
    Y_err[i] = np.sqrt(np.diag(D[i]))[1]
y=unp.uarray(Y,Y_err)

#Berechne Inhalte der Peaks mit der Integralformel fuer gauss
I=np.zeros(len(Y))
a=np.zeros(len(Y))
a_err=np.zeros(len(Y))
amp=np.zeros(len(Y))
amp_err=np.zeros(len(Y))
for i in range(len(Y)):
    a[i]=C[i][0]
    a_err[i]=np.sqrt(np.diag(D[i]))[0]
    amp[i]=C[i][3]
    amp_err[i]=np.sqrt(np.diag(D[i]))[3]
uamp=unp.uarray(amp,amp_err)

ua=unp.uarray(a,a_err)

I=uamp*unp.sqrt(np.pi)/unp.sqrt(a)
Zaehlrate=I/7764

#Berechne die Effizienz

omega=0.2

t=unp.uarray( 6084 ,1)

Aktivitaet=unp.uarray(4130,60)*unp.exp(-unp.log(2)/(unp.uarray(4943,5)*86400)*t)

W=np.array([28.6,7.6,0.4,26.5,2.2,3.1,2.0,0.9,12.9,4.2,14.6,0.6,10.2,13.6,1.6,21.0,0.5])/100

Q=4*np.pi/omega/Aktivitaet*Zaehlrate/W

print('Effizienz ',Q)

np.savetxt('Guess.txt',f(A,*params),delimiter=' ; ',newline=' ; ',fmt='%3d')

#Plots
#Potenzfunktion der Effizienz
def eff(x,p1,p2):
    return p1*np.power(x,p2)

Y_umgerechnet=y/params[0]-params[1]/params[0]

params3,covariance3 = curve_fit(eff,unp.nominal_values(Y_umgerechnet),unp.nominal_values(Q))
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 16

#plt.plot(X,unp.nominal_values(Q),'rx')
plt.plot(np.linspace(100,1500,100),eff(np.linspace(100,1500,100),*params3),label='Fitfunktion')
plt.errorbar(unp.nominal_values(Y_umgerechnet),unp.nominal_values(Q),xerr=unp.std_devs(Y_umgerechnet),yerr=unp.std_devs(Q),fmt='r.',label='Messergebnisse')
plt.ylabel('Effizienz')
plt.xlabel('Energie in eV')
plt.legend()
#plt.show()