import numpy as np
import matplotlib.pyplot as plt
import arquivo_heat as arqui
import matplotlib.ticker as ticker

def pos_thermal(path, n):
	
	k = 65
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
	k = 65
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
	dx = 0.02
	Tinf0 = 0.1
	Rex = np.zeros(n)
	Stx = np.zeros(n)
	for i in range(n):
		x = (x0 + i*dx)*L
		Rex[i] = (U*x)/nu
		Stx[i] = (dthmdy_m[i,0]*x/(Tinf0))/(Pr*Rex[i])
	
	return Rex, Stx

lambda_z = 'lam36'

path1 = '../SB/dx202/steady/%s/teste/' %(lambda_z)
path2 = '../SB2/dx202/unsteady/%s/fgv03/tt103680/' %(lambda_z)
path3 = '../SB2/dx202/unsteady/%s/fgv06/tt155520/' %(lambda_z)
path4 = '../SB2/dx202/unsteady/%s/fgv09/tt155520/' %(lambda_z)
path5 = '../SB2/dx202/unsteady/%s/fgv12/tt155520/' %(lambda_z)
path6 = '../SB2/dx202/unsteady/%s/fgv15/tt155520/' %(lambda_z)

full_path1 = path1 + 'heat_coeffs_00.dat'
print full_path1
data1 = open(full_path1)
x1, z1, dudy1, dthdy1, dubdy1, dthbdy1, dumdy1, dthmdy1 = arqui.arquivo_dat(data1)

x2, dthmdy2 = pos_thermal(path2, 1225)

x3, dthmdy3 = pos_thermal(path3, 1225)

x4, dthmdy4 = pos_thermal(path4, 1225)

x5, dthmdy5 = pos_thermal(path5, 1225)

x6, dthmdy6 = pos_thermal(path6, 1225)

k = 65
n1 = len(x1)/k 
n2 = len(x2)/k  
n3 = len(x3)/k 
n4 = len(x4)/k  
n5 = len(x5)/k  
n6 = len(x6)/k 

dthmdy1_m = matriz(dthmdy1, n1)
dthmdy2_m = matriz(dthmdy2, n2)
dthmdy3_m = matriz(dthmdy3, n3)
dthmdy4_m = matriz(dthmdy4, n4)
dthmdy5_m = matriz(dthmdy5, n5)
dthmdy6_m = matriz(dthmdy6, n6)

Rex1, Stx1 = stanton(dthmdy1_m, n1)
Rex2, Stx2 = stanton(dthmdy2_m, n2)
Rex3, Stx3 = stanton(dthmdy3_m, n3)
Rex4, Stx4 = stanton(dthmdy4_m, n4)
Rex5, Stx5 = stanton(dthmdy5_m, n5)
Rex6, Stx6 = stanton(dthmdy6_m, n6)

Pr = 0.72
x0 = 1.0
L = 0.1
dx = 0.02
U = 5.0
nu = 1.5094795344339218e-005
n7 = np.amax([n1, n2, n3, n4, n5, n6])
#n7 = n1
Rex7 = np.zeros(n7)
St_lam = np.zeros(n7)
St_tur = np.zeros(n7)

for i in range(n7):
	x = (x0 + i*dx)*L
	Rex7[i] = (U*x)/nu
	St_lam[i] = 0.33*(Pr**(-2.0/3.0))*(Rex7[i]**(-1.0/2.0))
	St_tur[i] = 0.0287*(Pr**(-2.0/3.0))*(Rex7[i]**(-1.0/5.0))

m = 1225
line = ['.','o', 'v', 's', '*', '^', 'p', 'd']
Rex_in = 0.0*L*U/nu
Rex_fi = 24.0*L*U/nu
markers_on = np.linspace(0,m-10,30,dtype='int')

@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%.2E" % x

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()

plt.plot(Rex7[0:n7], St_lam[0:n7],'--k', label='laminar')
plt.plot(Rex7[0:n7], St_tur[0:n7],'-.k', label='turbulento')
ax2.plot(Rex1[0:n1], Stx1[0:n1], marker='.', markevery=markers_on, label='GV_fgv00')
ax2.plot(Rex2[0:n2], Stx2[0:n2], marker='o', markevery=markers_on, label='GV_fgv03')
ax2.plot(Rex3[0:n3], Stx3[0:n3], marker='s', markevery=markers_on, label='GV_fgv06')
ax2.plot(Rex4[0:n4], Stx4[0:n4], marker='*', markevery=markers_on, label='GV_fgv09')
ax2.plot(Rex5[0:n5], Stx5[0:n5], marker='v', markevery=markers_on, label='GV_fgv12')
ax2.plot(Rex6[0:n6], Stx6[0:n6], marker='p', markevery=markers_on, label='GV_fgv15')
ax.set_xlabel('x')
ax.set_ylim(0.0005,0.0045)
ax.set_xlim(0,24.0)
ax.set_xticks(range(0,25,2))
ax2.set_xlabel('rex')
ax.set_ylabel('st')
ax2.set_xlim(Rex_in, Rex_fi)
ax2.set_xticks(np.linspace(Rex_in, Rex_fi, 5))
ax2.legend(loc=1, ncol=2, fontsize=8)
ax2.xaxis.set_major_formatter(major_formatter)
ax2.yaxis.set_major_formatter(major_formatter)

plt.show()

