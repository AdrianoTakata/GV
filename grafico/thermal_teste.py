import numpy as np
import matplotlib.pyplot as plt
import arquivo_heat as arqui

def pos_thermal(path, n):
	
	k = 33
	DthmDy = np.zeros(n*k)
	nome = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15']
	for i in range(16):
		name = 'heat_coeffs_%s.dat' %(nome[i])
		full_path = path + name
		data = open(full_path)
		x, z, dudy, dthdy, dubdy, dthbdy, dumdy, dthmdy = arqui.arquivo_dat(data)
		DthmDy = DthmDy + dthmdy
		
	return x, DthmDy/16.0

def matriz(dthmdy, n):

	cont = 0
	k = 33
	dthmdy_m = np.zeros((n,k))
	for j in range(k):
		for i in range(n):
			dthmdy_m[i,j] = dthmdy[cont]
			cont = cont +1
	return dthmdy_m

def stanton(dthmdy_m, n):
	
	Pr = 0.72
	U = 5.0
	nu = 1.5094795344339218e-005
	x0 = 1.0
	L = 0.1
	dx = 0.01
	Tinf0 = 0.1
	Rex = np.zeros(n)
	Stx = np.zeros(n)
	for i in range(n):
		x = (x0 + i*dx)*L
		Rex[i] = (U*x)/nu
		Stx[i] = (dthmdy_m[i,0]*x/(Tinf0))/(Pr*Rex[i])
	
	return Rex, Stx

lambda_z = 'lam18'

path1 = '../SB/steady/%s/teste1/' %(lambda_z)
path2 = '../SB/unsteady/%s/nssc/fgv03/tt172800/' %(lambda_z)
#path3 = '../SB/unsteady/%s/nssc/fgv06/tt172800/' %(lambda_z)
#path4 = '../SB/unsteady/%s/nssc/fgv09/tt172800/' %(lambda_z)
#path5 = '../SB/unsteady/%s/nssc/fgv12/tt120960/' %(lambda_z)
#path6 = '../SB/unsteady/%s/nssc/fgv15/tt120960/' %(lambda_z)

full_path1 = path1 + 'heat_coeffs.dat'
print full_path1
data1 = open(full_path1)
x1, z1, dudy1, dthdy1, dubdy1, dthbdy1, dumdy1, dthmdy1 = arqui.arquivo_dat(data1)

x2, dthmdy2 = pos_thermal(path2, 1865)

#x3, dthmdy3 = pos_thermal(path3, 1865)

#x4, dthmdy4 = pos_thermal(path4, 1865)

#x5, dthmdy5 = pos_thermal(path5, 1865)

#x6, dthmdy6 = pos_thermal(path6, 1865)

k = 33
n1 = len(x1)/k 
n2 = len(x2)/k  
#n3 = len(x3)/k 
#n4 = len(x4)/k  
#n5 = len(x5)/k  
#n6 = len(x6)/k 

dthmdy1_m = matriz(dthmdy1, n1)
dthmdy2_m = matriz(dthmdy2, n2)
#dthmdy3_m = matriz(dthmdy3, n3)
#dthmdy4_m = matriz(dthmdy4, n4)
#dthmdy5_m = matriz(dthmdy5, n5)
#dthmdy6_m = matriz(dthmdy6, n6)

Rex1, Stx1 = stanton(dthmdy1_m, n1)
Rex2, Stx2 = stanton(dthmdy2_m, n2)
#Rex3, Stx3 = stanton(dthmdy3_m, n3)
#Rex4, Stx4 = stanton(dthmdy4_m, n4)
#Rex5, Stx5 = stanton(dthmdy5_m, n5)
#Rex6, Stx6 = stanton(dthmdy6_m, n6)


Pr = 0.72
x0 = 1.0
L = 0.1
dx = 0.01
U = 5.0
nu = 1.5094795344339218e-005
#n7 = np.amax([n1, n2, n3, n4, n5, n6])
n7 = n2
Rex7 = np.zeros(n7)
St_lam = np.zeros(n7)
St_tur = np.zeros(n7)

for i in range(n7):
	x = (x0 + i*dx)*L
	Rex7[i] = (U*x)/nu
	St_lam[i] = 0.33*(Pr**(-2.0/3.0))*(Rex7[i]**(-1.0/2.0))
	St_tur[i] = 0.0287*(Pr**(-0.4))*(Rex7[i]**(-0.2))

m = n7
plt.figure()
plt.plot(Rex7[0:n7], St_lam[0:n7],'--k', label='laminar')
plt.plot(Rex7[0:n7], St_tur[0:n7],'-.k', label='turbulento')
plt.plot(Rex1[0:m], Stx1[0:m],'r', label='GV_fgv00')
plt.plot(X_tra, Y_rbf,'c', label='GV_fgv03')
#plt.plot(Rex3[0:m], Stx3[0:m],'m', label='GV_fgv06')
#plt.plot(Rex4[0:m], Stx4[0:m],'y', label='GV_fgv09')
#plt.plot(Rex5[0:n7], Stx5[0:n7],'g', label='GV_fgv12')
#plt.plot(Rex6[0:n7], Stx6[0:n7],'b', label='GV_fgv15')
plt.legend(loc=1, ncol=2)
#plt.ylim(0.0,0.0055)
plt.show()

