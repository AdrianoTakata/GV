import numpy as np
import matplotlib.pyplot as plt
import arquivo as arqui

Re = 33124.0

tipo = 'unsteady'
lam = 'lam09'

path1 = '../SB/%s/%s/nssc/fgv03/tt345600/' %(tipo, lam)
path2 = '../SB/%s/%s/nssc/fgv06/tt345600/' %(tipo, lam)
path3 = '../SB/%s/%s/nssc/fgv09/tt345600/' %(tipo, lam)
path4 = '../SB/%s/%s/nssc/fgv12/tt345600/' %(tipo, lam)
path5 = '../SB/%s/%s/nssc/fgv15/tt345600/' %(tipo, lam)

arq = 'mode00_02'

full_path1 = path1 + arq + '.dat'
print  full_path1 
data1 = open(full_path1)
x1, U_max1, alfa_ux1 = arqui.arquivo_dat(data1)
x1 = list(map(float, x1))
U_max1 = list(map(float, U_max1))
alfa_ux1 = list(map(float, alfa_ux1))
U_max1 = np.log10(U_max1)

full_path2 = path2 + arq + '.dat'
print full_path2
data2 = open(full_path2)
x2, U_max2, alfa_ux2 = arqui.arquivo_dat(data2)
x2 = list(map(float, x2))
U_max2 = list(map(float, U_max2))
alfa_ux2 = list(map(float, alfa_ux2))
U_max2 = np.log10(U_max2)

full_path3 = path3 + arq + '.dat'
print  full_path3
data3 = open(full_path3)
x3, U_max3, alfa_ux3 = arqui.arquivo_dat(data3)
x3 = list(map(float, x3))
U_max3 = list(map(float, U_max3))
alfa_ux3 = list(map(float, alfa_ux3))
U_max3 = np.log10(U_max3)

full_path4 = path4 + arq + '.dat'
print full_path4
data4 = open(full_path4)
x4, U_max4, alfa_ux4 = arqui.arquivo_dat(data4)
x4 = list(map(float, x4))
U_max4 = list(map(float, U_max4))
alfa_ux4 = list(map(float, alfa_ux4))
U_max4 = np.log10(U_max4)

full_path5 = path5 + arq + '.dat'
print full_path5
data5 = open(full_path5)
x5, U_max5, alfa_ux5 = arqui.arquivo_dat(data5)
x5 = list(map(float, x5))
U_max5 = list(map(float, U_max5))
alfa_ux5 = list(map(float, alfa_ux5))
U_max5 = np.log10(U_max5)
                        
plt.plot(x1, U_max1, label='fgv=03')
plt.plot(x2, U_max2, label='fgv=06')
plt.plot(x3, U_max3, label='fgv=09')
plt.plot(x4, U_max4, label='fgv=12')
plt.plot(x5, U_max5, label='fgv=15')
plt.xlabel('x')
plt.ylabel('U_max')
plt.xlim(0.0,14.0)
plt.ylim(-5.0,0.0)
plt.title(arq)
plt.legend(loc='best')
plt.show()
