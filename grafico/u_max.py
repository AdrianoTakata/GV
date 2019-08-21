import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui

arq = ['maxen00', 'maxen01', 'maxen02', 'maxen03', 'maxen04', 'maxen05', 'maxen06']
name = ['mode00', 'mode01', 'mode02', 'mode03', 'mode04', 'mode05', 'mode06']
path = '../steady/my_form=2/lambda_z=0036/imax=2825/'


for k in range(len(arq)):
	
	print 'arq = ', arq[k]
	full_path = path + arq[k] +'.dat'
	data = open(full_path)
	x, en, der = arqui.arquivo_dat(data)
	x = list(map(float, x))
	en = list(map(float, en))
	der = list(map(float, der))
        en = np.log10(en)

	plt.plot(x, en, label=name[k])
        plt.xlabel('x')
        plt.ylabel('log(Enk)')
        plt.xlim(4,11.5)
        plt.ylim(-10,-0)
	plt.legend(loc='best')

plt.show()

