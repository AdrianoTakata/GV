import sys
import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui
import matplotlib.ticker as ticker

arq = ['mode00_00', 'mode01_01', 'mode02_02', 'mode03_03', 'mode00_02', 'mode02_00', 'mode01_03', 'mode03_01',
      #'mode00_01', 'mode00_03', 'mode00_04', 'mode00_05', 'mode00_06','mode00_07', 'mode00_08', 'mode00_09', 'mode00_10'
      #'mode01_00', 'mode01_02', 'mode01_04', 'mode01_05', 'mode01_06', 'mode01_07', 'mode01_08', 'mode01_09', 'mode01_10'
      #'mode02_01', 'mode02_03', 'mode02_04', 'mode02_05', 'mode02_06', 'mode02_07', 'mode02_08', 'mode02_09', 'mode02_10'
      #'mode03_00', 'mode03_02', 'mode03_04', 'mode03_05', 'mode03_06', 'mode03_07', 'mode03_08', 'mode03_09', 'mode03_10'
      #'mode04_00', 'mode04_01', 'mode04_02', 'mode04_03', 'mode04_04', 'mode04_05', 'mode04_06', 'mode04_07', 'mode04_08', 'mode04_09', 'mode04_10'
      #'mode05_00', 'mode05_01', 'mode05_02', 'mode05_03', 'mode05_04', 'mode05_05', 'mode05_06', 'mode05_07', 'mode05_08', 'mode05_09', 'mode05_10'
      #'mode06_00', 'mode06_01', 'mode06_02', 'mode06_03', 'mode06_04', 'mode06_05', 'mode06_06', 'mode06_07', 'mode06_08', 'mode06_09', 'mode06_10'
      #'mode07_00', 'mode07_01', 'mode07_02', 'mode07_03', 'mode07_04', 'mode07_05', 'mode07_06', 'mode07_07', 'mode07_08', 'mode07_09', 'mode07_10'
      ]
#arq = ['mode00_00', 'mode01_01', 'mode00_02']
#arq = ['mode00_00','mode00_01','mode01_00']

name_arq = ['mode(0,0)', 'mode(1,1)', 'mode(2,2)', 'mode(3,3)', 'mode(0,2)', 'mode(2,0)','mode(1,3)','mode(3,1)',
           #'mode(0,1)', 'mode(0,3)', 'mode(0,4)', 'mode(0,5)', 'mode(0,6)', 'mode(0,7)', 'mode(0,8)', 'mode(0,9)', 'mode(0,10)'
           #'mode(1,0)', 'mode(1,2)', 'mode(1,4)', 'mode(1,5)', 'mode(1,6)', 'mode(1,7)', 'mode(1,8)', 'mode(1,9)', 'mode(1,10)'
           #'mode(2,1)', 'mode(2,3)', 'mode(2,4)', 'mode(2,5)', 'mode(2,6)', 'mode(2,7)', 'mode(2,8)', 'mode(2,9)', 'mode(2,10)'
           #'mode(3,0)', 'mode(3,2)', 'mode(3,4)', 'mode(3,5)', 'mode(3,6)', 'mode(3,7)', 'mode(3,8)', 'mode(3,9)', 'mode(3,10)'
           #'mode(4,0)', 'mode(4,1)', 'mode(4,2)', 'mode(4,3)', 'mode(4,4)', 'mode(4,5)', 'mode(4,6)', 'mode(4,7)', 'mode(4,8)', 'mode(4,9)', 'mode(4,10)'
           #'mode(5,0)', 'mode(5,1)', 'mode(5,2)', 'mode(5,3)', 'mode(5,4)', 'mode(5,5)', 'mode(5,6)', 'mode(5,7)', 'mode(5,8)', 'mode(5,9)', 'mode(5,10)'
           #'mode(6,0)', 'mode(6,1)', 'mode(6,2)', 'mode(6,3)', 'mode(6,4)', 'mode(6,5)', 'mode(6,6)', 'mode(6,7)', 'mode(6,8)', 'mode(6,9)', 'mode(6,10)'
           #'mode(7,0)', 'mode(7,1)', 'mode(7,2)', 'mode(7,3)', 'mode(7,4)', 'mode(7,5)', 'mode(7,6)', 'mode(7,7)', 'mode(7,8)', 'mode(7,9)', 'mode(7,10)'
           ]

line = ['.','o', 'v', 's', '*', '^', 'p', 'd']

Re = 33124.0
U = 5.0
nu = 1.5094795344339218e-005
x0 = 1.0
L = 0.1
dx = 0.02

lam = 'lam18'
base = 'nssc'
fgv = 'fgv06'
tt = 'tt207360'
dx = 'dx202'
ti = 'SB'
x_end = 20
y_end = 10

@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%.2E" % x

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
for k in range(len(arq)):

	path = '../%s/%s/unsteady/%s/%s/%s/' %(ti,dx, lam, fgv, tt)
#	path = '../SB/unsteady/lam18/teste/'
	full_path = path + arq[k] + '.dat'
	print 'arq =', arq
	data = open(full_path)
	x, U_max, alfa_ux = arqui.arquivo_dat(data)
	x = list(map(float, x))
	U_max = list(map(float, U_max))
	alfa_ux = list(map(float, alfa_ux))
        U_max = np.log10(U_max)
	print len(x)
        markers_on = np.linspace(0, len(x)-4, 25, dtype='int', endpoint=True)
	NN = len(x)
	Rex1 = 4.0*L*U/nu
	Rex2 = x_end*L*U/nu
		
	if (k < 8):     
        	ax.plot(x, U_max, marker=line[k], markevery=markers_on, label=name_arq[k])
#		ax.plot(x, U_max)
	else:
		ax.plot(x, U_max, marker=line[k], markevery=markers_on, label=name_arq[k])
#		ax.plot(x, U_max)

        ax.set_xlabel('x')
        ax.set_ylabel('log(Umax)')
        ax.set_xlim(4.0,x_end)
	ax.set_ylim(-y_end,0.0)
	ax.set_xticks(range(4,x_end+1,2))
	ax.set_yticks(range(-y_end, 1,2))
	ax2.set_xlabel('Rex')
	ax2.set_xlim(Rex1, Rex2)
	ax2.set_xticks(np.linspace(Rex1, Rex2, 7))
	ax2.xaxis.set_major_formatter(major_formatter)
	ax.legend(loc=4, ncol=2, fontsize=8)#, mode='expand')

plt.show()
print Rex1
print Rex2

