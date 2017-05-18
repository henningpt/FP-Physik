import numpy as np 
f600,f525,f450,f375,f300,f225,f150,f75=np.genfromtxt('4_1_peaks.txt', unpack=True)

D600=np.zeros(shape=(8,1))
for i in range(7):
    D600[i]=-(f600[i]-f600[i+1])

D525=np.zeros(shape=(8,1))
for i in range(7):
    D525[i]=-(f525[i]-f525[i+1])

D450=np.zeros(shape=(8,1))
for i in range(7):
    D450[i]=-(f450[i]-f450[i+1])

D375=np.zeros(shape=(8,1))
for i in range(7):
    D375[i]=-(f375[i]-f375[i+1])

D300=np.zeros(shape=(8,1))
for i in range(7):
    D300[i]=-(f300[i]-f300[i+1])
    
D225=np.zeros(shape=(8,1))
for i in range(7):
    D225[i]=-(f225[i]-f225[i+1])

D150=np.zeros(shape=(8,1))
for i in range(7):
    D150[i]=-(f150[i]-f150[i+1])

D75=np.zeros(shape=(8,1))
for i in range(7):
    D75[i]=-(f75[i]-f75[i+1])

print(D300)
print([D300.T,D600.T])

np.savetxt('4_1_diff.txt', [D600,D525].T,delimiter='&')
