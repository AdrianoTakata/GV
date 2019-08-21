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
epsilon = 5.0e-4
Re_lambda = U*lam/nu
rt = epsilon*Re_lambda
#rt = 0.01
print rt

arq = [#'mode00_00', 'mode01_01', 'mode02_02', 'mode03_03']
        'mode02_00', 'mode00_02', 'mode01_03', 'mode03_01']

name_arq = [#'mode(0,0)', 'mode(1,1)', 'mode(2,2)', 'mode(3,3)'
            'mode(0,2)', 'mode(2,0)','mode(1,3)','mode(3,1)']

path1 = '../dados/resul_marensi/'
path2 = '../marensi/dx202/A504/'

for k in range(len(arq)):
	
	full_path1 = path1 + arq[k] + '.dat'
	print full_path1
	data1  = np.loadtxt(full_path1, dtype='str')
	x1 = data1[:,0]
	en1 = data1[:,1]
	for i in range(len(x1)):
		x1[i]  = x1[i].replace(',','.')
        	en1[i] = en1[i].replace(',','.')

	x1 = list(map(float, x1))
	en1 = list(map(float, en1))
	for i in range(len(x1)):
		x1[i] = x1[i]

	
	full_path2 = path2 + arq[k] + '.dat'
	print full_path2
	#full_path2 = path2 + 'en01.dat'
	data2 = open(full_path2)
	x, U_max, alfa_ux = arqui.arquivo_dat(data2)
	x = list(map(np.float64, x))
	U_max = list(map(np.float64, U_max)) 
	alfa_ux = list(map(np.float64, alfa_ux))

	for t in range(len(U_max)):
		if k == 0 or k == 4:
			U_max[t] = (2.0*rt*U_max[t]**2.0)
		else:
			U_max[t] = ((1.0/2.0)*rt*U_max[t]**2.0)
	U_max = np.log10(U_max)

	for t in range(len(x)):
		#x[t] = x[t]**2.0/Re # x admensional
		x[t] = x[t]*L # x_star dimensional
		x[t] = x[t]/lam 
		x[t] = x[t]*kx
	
	if k == 0 or k == 4:
                plt.plot(x1, en1, '.b' )#,label='Marensi & Ricco, 2017')
                plt.plot(x, U_max, label=name_arq[k])
	else:
		plt.plot(x1, en1, '.b')
                plt.plot(x, U_max, label=name_arq[k])


plt.xlabel(r'$\bar x $')
plt.ylabel(r'$log(Umax)$')
plt.xlim(3.0,8.0)
plt.ylim(-16,0)
plt.yticks(range(-16,1,2))
plt.legend(loc='best', ncol=2, fontsize=8)
plt.show()

