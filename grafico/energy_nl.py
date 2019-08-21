import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui

arq = ['en00.dat', 'en01.dat', 'en02.dat', 'en03.dat', 'en04.dat', 'en05.dat', 'en06.dat']

lambda_z = 'lam18'

#path1 = '../../../Dropbox/Adriano_Takata/dados/'
path1 = '../dados/Stacionary/%s/' %(lambda_z)
path2 = '../SB/dx102/steady/%s/nssc/' %(lambda_z)
path3 = '../SB/dx202/steady/%s/nssc/' %(lambda_z)
path4 = '../SB/dx102/steady/%s/nssc/' %(lambda_z)

x_ini = 4.0
x_end = 11.5
y_ini = -16.0
y_end = 0.0

for k in range(len(arq)):
	
	print 'arq = ', arq[k]
	'''
	full_path1 = path1 + arq[k]
	print 'arq =', arq[k]
	data1 = open(full_path1)
	x1, en1, der1 = arqui.arquivo_dat(data1)
	x1 = list(map(float, x1))
	en1 = list(map(float, en1)) 
	der1 = list(map(float, der1))
	'''
	full_path1 = path1 + arq[k]
	data1  = np.loadtxt(full_path1, dtype='str')
	x1 = data1[:,0]
	en1 = data1[:,1]
	for i in range(len(x1)):
	
		x1[i]  = x1[i].replace(',','.')
		en1[i] = en1[i].replace(',','.')

	x1 = list(map(float, x1))
	en1 = list(map(float, en1)) 
        
	full_path2 = path2 + arq[k]
	data2 = open(full_path2)
	x2, en2, der2 = arqui.arquivo_dat(data2)
	x2 = list(map(float, x2))
	en2 = list(map(float, en2))
	der2 = list(map(float, der2))
        	
        full_path3 = path3 + arq[k]
	data3 = open(full_path3)
	x3, en3, der3 = arqui.arquivo_dat(data3)
	x3 = list(map(float, x3))
	en3 = list(map(float, en3))
	der = list(map(float, der3))
	
        full_path4 = path4 + arq[k]
	data4 = open(full_path4)
	x4, en4, der4 = arqui.arquivo_dat(data4)
	x4 = list(map(float, x4))
	en4 = list(map(float, en4))
	der4 = list(map(float, der4))
	
	
	if k == 0:
		plt.plot(x1, en1,'.b',label='Souza, 2017')
#		plt.plot(x2, en2,'--r', label='blasius')
		plt.plot(x3, en3,'-k', label='Simulation')
		plt.plot(x4, en4,'-g', label='nscc')
 			
	else:
		plt.plot(x1, en1, '.b')# ,label='Reference')
#		plt.plot(x2, en2,'--r')# , label='Obtained')
		plt.plot(x3, en3,'-k')
		plt.plot(x4, en4,'-g')

plt.xlabel('x')
plt.ylabel('log(Enk)')
plt.xlim(x_ini,x_end)
plt.ylim(y_ini,y_end)
plt.legend(loc='lower right')
plt.show()

