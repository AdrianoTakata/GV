import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

path1 = 'lam36_fgv15_perd00.dat'
teste1 = np.loadtxt(path1)

I = 2425
J = 257
K = 32

Y = np.zeros((J,K+1))
Z = np.zeros((J,K+1))
UX = np.zeros((J,K+1))
TH = np.zeros((J,K+1))
Vp = np.zeros((J,K+1))
Wp = np.zeros((J,K+1))

cont = 0
for k in range(K):
	for j in range(J):
		Y[j,k] = teste1[cont,0]
		Z[j,k] = teste1[cont,1]
		UX[j,k] = teste1[cont,2]
		TH[j,k] = teste1[cont,3]
		Vp[j,k] = teste1[cont,4]
		Wp[j,k] = teste1[cont,5]
		cont += 1
		if (cont >= I*J*K):
			break

Z[:,K] = 0.36*np.ones(J)
Y[:,K] = Y[:,0]
UX[:,K] = UX[:,0]
TH[:,K] = TH[:,0]
Vp[:,K] = Vp[:,0]
Wp[:,K] = Wp[:,0]

levels = np.linspace(0.05, 0.95, 10)

Z1 = Z[0:J-1:4,:]
Y1 = Y[0:J-1:4,:]
UX1 = UX[0:J-1:4,:]
TH1 = TH[0:J-1:4,:]
Vp1 = Vp[0:J-1:4,:]
Wp1 = Wp[0:J-1:4,:]

plt.figure()
CS1 = plt.contour(Z1, Y1, UX1, levels, cmap=plt.cm.jet)
CS2 = plt.contour(Z1, Y1, TH1, levels, cmap=plt.cm.jet, linestyles='dashed')
q = plt.quiver(Z1, Y1, Wp1, Vp1)
plt.ylim(0.0, 0.24)
plt.xlabel('z')
plt.ylabel('y')
#plt.xlim(0.0, 0.09)
#plt.clabel(CS1, inline=1)
#plt.clabel(CS2, inline=1)
plt.show()

