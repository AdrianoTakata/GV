import sys
import numpy as np
import matplotlib.pyplot as plt
import arquivo_sec as arqui
import matplotlib.ticker as ticker

arq = ['mode0_0', 'mode00',
       'mode01', 'mode02', 'mode03', 'mode04', 'mode05', 'mode06',
       'mode07', 'mode08', 'mode09', 'mode10', 'mode11', 'mode12', 'mode13','mode14','mode15','mode16'
      ]

name_arq = ['mode0_0', 'mode00',
       'mode01', 'mode02', 'mode03', 'mode04', 'mode05', 'mode06',
       'mode07', 'mode08', 'mode09', 'mode10', 'mode11', 'mode12', 'mode13','mode14','mode15','mode16'
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
tt = 'sino'
dx = 'dx202'
ti = 'SB'
x_end = 20
y_end = 6

@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%.2E" % x

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
for k in range(len(arq)):

	path = '../%s/%s/sec/%s/%s/' %(ti,dx, lam, tt)
#	path = '../SB/unsteady/lam18/teste/'
	full_path = path + arq[k] + '.dat'
	print 'arq =', arq
	data = open(full_path)
	x, U_max, alfa_x, U_maxd2, alpha_xd2 = arqui.arquivo_dat(data)
	x = list(map(float, x))
	U_max = list(map(float, U_max))
	alfa_x = list(map(float, alfa_x))
	U_maxd2 = list(map(float, U_maxd2))
	alpha_xd2 = list(map(float, alpha_xd2))
        U_max = np.log10(U_max)
	U_maxd2 = np.log10(U_maxd2)
	print len(x)
        markers_on = np.linspace(0, len(x)-4, 25, dtype='int', endpoint=True)
	NN = len(x)
	Rex1 = 4.0*L*U/nu
	Rex2 = x_end*L*U/nu
		
	if (k < len(arq)):     
#        	ax.plot(x, U_max, marker=line[k], markevery=markers_on, label=name_arq[k])
#		ax.plot(x, U_max)
		ax.plot(x, U_maxd2)
	else:
#		ax.plot(x, U_max, marker=line[k], markevery=markers_on, label=name_arq[k])
#		ax.plot(x, U_max)
		ax.plot(x, U_maxd2)

        ax.set_xlabel('x')
        ax.set_ylabel('log(Umaxd2)')
        ax.set_xlim(4.0,x_end)
	ax.set_ylim(-y_end,0.0)
	ax.set_xticks(range(4,x_end+1,2))
	ax.set_yticks(range(-y_end, 1,2))
	ax2.set_xlabel('Rex')
	ax2.set_xlim(Rex1, Rex2)
#	ax2.set_xticks(np.linspace(Rex1, Rex2, 7))
	ax2.xaxis.set_major_formatter(major_formatter)
	ax.legend(loc=4, ncol=2, fontsize=8)#, mode='expand')

plt.show()
print Rex1
print Rex2

