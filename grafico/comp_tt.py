import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui
 
arq = ['mode00_00', 'mode01_01', 'mode02_02', 'mode03_03', 'mode00_02', 'mode02_00','mode01_03','mode03_01'] 
#arq = ['mode00_00', 'mode01_01', 'mode02_02', 'mode03_03']
#arq = ['mode00_00', 'mode01_01', 'mode00_02']

Re = 33124.0

lam = 'lam18'
fgv = 'fgv12'
t1 = 'tt103680'
t2 = 'tt155520'

for k in range(len(arq)):

	path = '../SB/dx202/unsteady/%s/%s/%s/' %(lam, fgv, t1)
	full_path = path + arq[k] + '.dat'
	print 'arq =', arq
	data = open(full_path)
	x, U_max, alfa_ux = arqui.arquivo_dat(data)
	x = list(map(float, x))
	U_max = list(map(float, U_max))
	alfa_ux = list(map(float, alfa_ux))
 
        U_max = np.log10(U_max)
       
	path1 = '../SB/dx202/unsteady/%s/%s/%s/' %(lam, fgv, t2)
	full_path1 = path1 + arq[k] + '.dat'
	print 'arq =', arq
	data1 = open(full_path1)
	x1, U_max1, alfa_ux1 = arqui.arquivo_dat(data1)
	x1 = list(map(float, x1))
	U_max1 = list(map(float, U_max1))
	alfa_ux1 = list(map(float, alfa_ux1))
 
        U_max1 = np.log10(U_max1)
 	        
        plt.plot(x, U_max, 'k', label='SB')
	plt.plot(x1, U_max1,'r' ,label='SB2')
        plt.xlabel('x')
        plt.ylabel('log(U_max)')
        plt.xlim(0.0,25.0)
        plt.ylim(-12,-0)
        plt.legend(loc='best')

plt.show()

