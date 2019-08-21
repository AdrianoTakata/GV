import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui
 
Re = 60240.039132584265
U = 9.18
L = 0.1
nu = 1.523903392525267e-05
lambda_z = 0.008
lam = lambda_z/(2.0*np.pi)
kx = 4.35e-3
epsilon = 5.3e-4
Re_lambda = U*lam/nu
rt = epsilon*Re_lambda
#rt = 0.01
print rt


arq = ['mode00_00', 'mode01_01','mode02_02', 'mode03_03','mode02_00', 'mode00_02','mode01_03','mode03_01']

path1 = '../marensi/unsteady/lam08/teste/tt92000/'
path2 = '../marensi/unsteady/lam08/teste/tt184000/'

for k in range(len(arq)):

	full_path1 = path1 + arq[k] + '.dat'
	print full_path1
	data1 = open(full_path1)
	x1, U_max1, alfa_ux1 = arqui.arquivo_dat(data1)
	x1 = list(map(np.float64, x1))
	U_max1 = list(map(np.float64, U_max1)) 
	alfa_ux1 = list(map(np.float64, alfa_ux1))

	for t in range(len(U_max1)):
		if k == 0 and k == 4:
			U_max1[t] = (2.0*rt*U_max1[t]**2.0)
		else:
			U_max1[t] = ((1.0/2.0)*rt*U_max1[t]**2.0)

	U_max1 = np.log10(U_max1)

	for t in range(len(x1)):
		x1[t] = x1[t]*L # x_star dimensional
		x1[t] = x1[t]/lam 
		x1[t] = x1[t]*kx
	

		
	full_path2 = path2 + arq[k] + '.dat'
	print full_path2
	data2 = open(full_path2)
	x2, U_max2, alfa_ux2 = arqui.arquivo_dat(data2)
	x2 = list(map(np.float64, x2))
	U_max2 = list(map(np.float64, U_max2)) 
	alfa_ux2 = list(map(np.float64, alfa_ux2))

	for t in range(len(U_max2)):
		if k == 0 and k==4 :
			U_max2[t] = (2.0*rt*U_max2[t]**2.0)
		else:
			U_max2[t] = ((1.0/2.0)*rt*U_max2[t]**2.0)

	U_max2 = np.log10(U_max2)

	for t in range(len(x2)):
		x2[t] = x2[t]*L # x_star dimensional
		x2[t] = x2[t]/lam 
		x2[t] = x2[t]*kx
	
	if k == 0 and k == 4:
                plt.plot(x1, U_max1, '--r' , label='1')
                plt.plot(x2, U_max2,'-.k', label='2')
#                plt.xlabel(r'$\bar x$',fontsize=14)
#                plt.ylabel('log(Umax)',fontsize=12)
#                plt.xlim(0.0,12.0)
#                plt.ylim(-16.0,-0)
#                plt.legend(loc='best')
        else:
                plt.plot(x1, U_max1, '--r',label='1')
                plt.plot(x2, U_max2,'-.k', label='2')
                plt.xlabel(r'$\bar x$',fontsize=14)
                plt.ylabel('log(Umax)',fontsize=12)
		plt.xlim(3.0,8.0)
                plt.ylim(-16,-0)
                plt.legend(loc='best')

plt.show()

