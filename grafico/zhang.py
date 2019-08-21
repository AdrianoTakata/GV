import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui

Re = 31250.0
epsilon1 = 1.0e-3
epsilon2 = 7.0e-4
epsilon3 = 7.0e-4
epsilon4 = 6.0e-4
nu = 1.6e-5
lambda_z = 0.023
U = 5.0
Re_lambda = (U*lambda_z)/(2.0*np.pi*nu)
#Re_lambda = 1145.0
print Re_lambda
rt1 = epsilon1*Re_lambda
rt2 = epsilon2*Re_lambda
rt3 = epsilon3*Re_lambda
rt4 = epsilon4*Re_lambda
L = 0.1


arq = ['mode00_00', 'mode01_01', 'mode02_02']#, 'mode02_00', 'mode00_02']
#arq = ['mode02_00', 'mode00_02'] 
name_arq = ['mode(0,0)', 'mode(1,1)', 'mode(2,2)', 'mode(2,0)', 'mode(0,2)']

path = '../dados/resul_zhang/f=10hz/'
path1 = '../zhang/fgv10/semdis2/'
path2 = '../zhang/fgv10/tt103680/'
#path3 = '../zhang/unsteady/fgv01/i1850/'
#path4 = '../zhang/unsteady/fgv01/i1749/'

for k in range(len(arq)):
	
	full_path = path + arq[k] + '.dat'
	print full_path
	data  = np.loadtxt(full_path, dtype='str')
	x = data[:,0]
	en = data[:,1]
	for i in range(len(x)):
		x[i]  = x[i].replace(',','.')
        	en[i] = en[i].replace(',','.')

	x = list(map(float, x))
	en = list(map(float, en))
	
	full_path1 = path1 + arq[k] + '.dat'
	print full_path1
	data1 = open(full_path1)
	x1, U_max1, alfa_ux1 = arqui.arquivo_dat(data1)
	x1 = list(map(np.float64, x1))
	U_max1 = list(map(np.float64, U_max1)) 
	alfa_ux1 = list(map(np.float64, alfa_ux1))
	
	full_path2 = path2 + arq[k] + '.dat'
	print full_path2
	data2 = open(full_path2)
	x2, U_max2, alfa_ux2 = arqui.arquivo_dat(data2)
	x2 = list(map(np.float64, x2))
	U_max2 = list(map(np.float64, U_max2)) 
	alfa_ux2 = list(map(np.float64, alfa_ux2))
	'''	
	full_path3 = path3 + arq[k] + '.dat'
	print full_path3
	data3 = open(full_path3)
	x3, U_max3, alfa_ux3 = arqui.arquivo_dat(data3)
	x3 = list(map(np.float64, x3))
	U_max3 = list(map(np.float64, U_max3)) 
	alfa_ux3 = list(map(np.float64, alfa_ux3))

	full_path4 = path4 + arq[k] + '.dat'
	print full_path4
	data4 = open(full_path4)
	x4, U_max4, alfa_ux4 = arqui.arquivo_dat(data4)
	x4 = list(map(np.float64, x4))
	U_max4 = list(map(np.float64, U_max4)) 
	alfa_ux4 = list(map(np.float64, alfa_ux4))
	'''
	x_hat = np.zeros((len(x1),1))
	for t in range(len(x1)):
		#x[t] = (x[t]**2.0)/Re
		x1[t] = x1[t]*L 
		x_hat[t] = (2.0*np.pi*x1[t]) / (Re_lambda*lambda_z)
		if k == 0 or k == 3:
			U_max1[t] = np.log10(1.05*rt1*U_max1[t])
			U_max2[t] = np.log10(1.6*rt2*U_max2[t])
#			U_max3[t] = np.log10(1.6*rt3*U_max3[t])
#			U_max4[t] = np.log10(1.8*rt4*U_max4[t])

		else:
			U_max1[t] = np.log10(0.525*rt1*U_max1[t])
			U_max2[t] = np.log10(0.775*rt2*U_max2[t])
#			U_max3[t] = np.log10(0.775*rt3*U_max3[t])
#			U_max4[t] = np.log10(0.875*rt4*U_max4[t])


	if k == 0:
                plt.plot(x, en, '.b' )#,label='Marensi & Ricco, 2017')
                plt.plot(x_hat, U_max1, '-k',label='normal')
		plt.plot(x_hat, U_max2, '-r', label='i1951')
#		plt.plot(x_hat, U_max3, '-g',label='i1850')
		#plt.plot(x_hat, U_max4, '-m',label='i1749')

	else:
		plt.plot(x, en, '.b')
                plt.plot(x_hat, U_max1, '-k')
		plt.plot(x_hat, U_max2, '-r')
#		plt.plot(x_hat, U_max3, '-g')
		#plt.plot(x_hat, U_max4, '-m')

plt.xlabel('x')
plt.ylabel('En')
plt.xlim(0.05,0.4)
plt.ylim(-4.0,-0.0)
plt.legend(loc='best')

plt.show()

