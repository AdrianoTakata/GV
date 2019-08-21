import numpy as np
import matplotlib.pyplot as plt
import arquivo_heat as arqui
import matplotlib.ticker as ticker


def pos_cfx(path, n):
	
	k = 65
	DumDy = np.zeros(n*k)
	nome = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15']
	for i in range(16):
		name = 'heat_coeffs_%s.dat' %(nome[i])
		full_path = path + name
		data = open(full_path)
		x, z, dudy, dthdy, dubdy, dthbdy, dumdy, dthmdy = arqui.arquivo_dat(data)
		DumDy = DumDy + dumdy
		
	return x, DumDy/16.0

def matriz(dumdy, n):

	cont = 0
	k = 65
	dumdy_m = np.zeros((n,k))
	for j in range(k):
		for i in range(n):
			dumdy_m[i,j] = dumdy[cont]
			cont = cont +1
	return dumdy_m

def cf(dumdy_m, n):
	
	rho = 1315.0
	U = 5.0
	nu = 1.5094795344339218e-005
	x0 = 1.0
	L = 0.1
	dx = 0.02
	Rex = np.zeros(n)
	Cfx = np.zeros(n)
	for i in range(n):
		x = (x0 + i*dx)*L
		Rex[i] = (U*x)/nu
		Cfx[i] = dumdy_m[i,0]/(0.5*rho*(U**2))
	
	return Rex, Cfx

lambda_z = 'lam18'

path1 = '../SB/dx202/steady/%s/nssc/' %(lambda_z)
path2 = '../SB/dx202/unsteady/%s/fgv03/tt155520/' %(lambda_z)
path3 = '../SB/dx202/unsteady/%s/fgv06/tt207360/' %(lambda_z)
path4 = '../SB/dx202/unsteady/%s/fgv09/tt207360/' %(lambda_z)
path5 = '../SB/dx202/unsteady/%s/fgv12/tt207360/' %(lambda_z)
path6 = '../SB/dx202/unsteady/%s/fgv15/tt207360/' %(lambda_z)

full_path1 = path1 + 'heat_coeffs_00.dat'
print full_path1
data1 = open(full_path1)
x1, z1, dudy1, dthdy1, dubdy1, dthbdy1, dumdy1, dthmdy1 = arqui.arquivo_dat(data1)

x2, dumdy2 = pos_cfx(path2, 1065)

x3, dumdy3 = pos_cfx(path3, 1065)

x4, dumdy4 = pos_cfx(path4, 1065)

x5, dumdy5 = pos_cfx(path5, 1065)

x6, dumdy6 = pos_cfx(path6, 1065)

k = 65
n1 = len(x1)/k 
n2 = len(x2)/k  
n3 = len(x3)/k 
n4 = len(x4)/k  
n5 = len(x5)/k  
n6 = len(x6)/k 

dumdy1_m = matriz(dumdy1, n1)
dumdy2_m = matriz(dumdy2, n2)
dumdy3_m = matriz(dumdy3, n3)
dumdy4_m = matriz(dumdy4, n4)
dumdy5_m = matriz(dumdy5, n5)
dumdy6_m = matriz(dumdy6, n6)

Rex1, Cfx1 = cf(dumdy1_m, n1)
Rex2, Cfx2 = cf(dumdy2_m, n2)
Rex3, Cfx3 = cf(dumdy3_m, n3)
Rex4, Cfx4 = cf(dumdy4_m, n4)
Rex5, Cfx5 = cf(dumdy5_m, n5)
Rex6, Cfx6 = cf(dumdy6_m, n6)

Pr = 0.72
x0 = 1.0
L = 0.1
dx = 0.02
U = 5.0
nu = 1.5094795344339218e-005
n7 = np.amax([n1, n2, n3, n4, n5])
#n7 = n4
Rex7 = np.zeros(n7)
Cfx_lam = np.zeros(n7)
Cfx_tur = np.zeros(n7)

for i in range(n7):
	x = (x0 + i*dx)*L
	Rex7[i] = (U*x)/nu
	Cfx_lam[i] = 0.664*(Rex7[i]**(-1.0/2.0))
	Cfx_tur[i] = 0.027*(Rex7[i]**(-1.0/7.0))

m = 1065
line = ['.','o', 'v', 's', '*', '^', 'p', 'd']
Rex_in = 0.0*L*U/nu
Rex_fi = 20.0*L*U/nu
markers_on = np.linspace(0, m-10, 30, dtype='int', endpoint=True)
'''
xp = 650

print 'St in x=14.0'
print 'fgv=0 -> %e' %(Stx1[xp])
print 'fgv=3 -> %e' %(Stx2[xp])
print 'fgv=6 -> %e' %(Stx3[xp])
print 'fgv=9 -> %e' %(Stx4[xp])
print 'fgv=12 -> %e' %(Stx5[xp])
print 'fgv=15 -> %e' %(Stx6[xp])

print 'laminar'
print 'fgv=0', (Stx1[xp]*100.0)/St_lam[xp] - 100.0
print 'fgv=3', (Stx2[xp]*100.0)/St_lam[xp] - 100.0
print 'fgv=6', (Stx3[xp]*100.0)/St_lam[xp] - 100.0
print 'fgv=9', (Stx4[xp]*100.0)/St_lam[xp] - 100.0
print 'fgv=12', (Stx5[xp]*100.0)/St_lam[xp] - 100.0
print 'fgv=15', (Stx6[xp]*100.0)/St_lam[xp] - 100.0

print 'turbulento'
print 'fgv=0', 100.0 - (Stx1[xp]*100.0)/St_tur[xp]
print 'fgv=3', 100.0 - (Stx2[xp]*100.0)/St_tur[xp] 
print 'fgv=6', 100.0 - (Stx3[xp]*100.0)/St_tur[xp] 
print 'fgv=9', 100.0 - (Stx4[xp]*100.0)/St_tur[xp] 
print 'fgv=12', 100.0 - (Stx5[xp]*100.0)/St_tur[xp]
print 'fgv=15', 100.0 - (Stx6[xp]*100.0)/St_tur[xp] 
'''
@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%.2E" % x

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
ax2.plot(Rex7[0:m], Cfx_lam[0:m],'--k', label='laminar')
ax2.plot(Rex7[0:m], Cfx_tur[0:m],'-.k', label='turbulent')
ax2.plot(Rex1[0:m], Cfx1[0:m], marker='.', markevery=markers_on, label='GV_fgv00')
ax2.plot(Rex2[0:m], Cfx2[0:m], marker='o', markevery=markers_on, label='GV_fgv03')
ax2.plot(Rex3[0:m], Cfx3[0:m], marker='s', markevery=markers_on, label='GV_fgv06')
ax2.plot(Rex4[0:m], Cfx4[0:m], marker='*', markevery=markers_on, label='GV_fgv09')
ax2.plot(Rex5[0:m], Cfx5[0:m], marker='v', markevery=markers_on, label='GV_fgv12')
ax2.plot(Rex6[0:m], Cfx6[0:m], marker='p', markevery=markers_on, label='GV_fgv15')
ax.set_xlabel('x')
ax.set_ylim(0.0005,0.0065)
ax.set_xlim(0,20.0)
ax.set_xticks(range(0,21,2))
#ax.set_yticks(np.linspace(0.0005,0.0055))
ax2.set_xlabel('Rex')
ax.set_ylabel('Cf')
ax2.set_xlim(Rex_in, Rex_fi)
ax2.set_xticks(np.linspace(Rex_in, Rex_fi, 7))
ax2.legend(loc=1, ncol=2, fontsize=8)
ax2.xaxis.set_major_formatter(major_formatter)
ax2.yaxis.set_major_formatter(major_formatter)

plt.show()

