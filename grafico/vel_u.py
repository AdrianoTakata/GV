import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

U = 5.0
nu = 1.5094795344339218e-005
x0 = 1.0
L = 0.1
pos_x = 1270
dx = 0.01
x = (x0 * pos_x*dx)*L
I = 1625
J = 257
K = 32

path1 = 'lam09_blasius.dat'
data1 = np.loadtxt(path1)
X = np.zeros((I,J))
Y = np.zeros((I,J))
UX = np.zeros((I,J))
TH = np.zeros((I,J))

cont = 0
for j in range(J):
	for i in range(I):
		X[i,j] = data1[cont,0]
		Y[i,j] = data1[cont,1]
		UX[i,j] = data1[cont,2]
		cont += 1
		if (cont >= I*J):
			break

y = Y[1,:]
y = y*L
eta = y*np.sqrt(U/(2.0*nu*x))
plt.plot(UX[pos_x,:], eta,label='blasius')

for i in np.array([0, 4, 8, 12]):
	if i == 12:
		path2 = 'lam09_fgv03_perd%s.dat' %(i)
	else:
		path2 = 'lam09_fgv03_perd0%s.dat' %(i)
	data2 = np.loadtxt(path2)


	Y = np.zeros((J,K))
	Z = np.zeros((J,K))
	UX = np.zeros((J,K))
	TH = np.zeros((J,K))

	cont = 0
	for k in range(K):
		for j in range(J):
			Y[j,k] = data2[cont,0]
			Z[j,k] = data2[cont,1]
			UX[j,k] = data2[cont,2]
			cont += 1
			if (cont >= I*J*K):
				break

	plt.plot(UX[:,15],eta,label=i)
	plt.ylim(0.0,7.0)
	plt.xlim(0.0,1.2)

plt.xlabel('U')
plt.ylabel(r'$\eta$')
plt.legend()
plt.show()
