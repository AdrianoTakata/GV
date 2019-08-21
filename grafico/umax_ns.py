import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui

path = '../unsteady/lambda_z=0009/omega_gv=03/imax=2825/'

full_path = path + 'ns_rms.dat'
data = open(full_path)
x, en, der = arqui.arquivo_dat(data)
x = list(map(float, x))
en = list(map(float, en))
der = list(map(float, der))
en = np.log10(en)


plt.plot(x, en,'k')
plt.xlabel('x')
plt.ylabel('log(Enk)')
#plt.xlim(4,16.0)
#plt.ylim(-8,-0)
plt.legend(loc='best')

plt.show()

