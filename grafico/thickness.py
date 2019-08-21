import numpy as np
import matplotlib.pyplot as plt
import arquivo3 as arqui3


U = 9.18
nu = 1.523903392525267e-05
L = 0.1
Re = (U*L)/nu

path1 = '../../../Dropbox/Doutorado/baseflow/'
path2 = '../marensi/steady/grad00/'
path3 = '../marensi/steady/nssc/'
path4 = '../marensi/steady/grad_nssc/'
path5 = '../marensi/steady/grad003/'

full_path1 = path1 + 'delta1.dat'
data1 = np.loadtxt(full_path1, dtype='str')
print data1
x1 = data1[:,0]
delta1 = data1[:,1]
delta2 = data1[:,2]
H12 = data1[:,3]

full_path2 = path2 + 'H12.dat'
data2 = open(full_path2)
x2, delta_2, delta1_2, delta2_2, H12_2 = arqui3.arquivo_dat(data2)
x2 = list(map(float, x2))
delta1_2 = list(map(float, delta1_2))
delta2_2 = list(map(float, delta2_2))
H12_2 = list(map(float, H12_2))

full_path3 = path3 + 'H12.dat'
data3 = open(full_path3)
x3, delta_3, delta1_3, delta2_3, H12_3 = arqui3.arquivo_dat(data3)
x3 = list(map(float, x3))
delta1_3 = list(map(float, delta1_3))
delta2_3 = list(map(float, delta2_3))
H12_3 = list(map(float, H12_3))

full_path4 = path4 + 'H12.dat'
data4 = open(full_path4)
x4, delta_4, delta1_4, delta2_4, H12_4 = arqui3.arquivo_dat(data4)
x4 = list(map(float, x4))
delta1_4 = list(map(float, delta1_4))
delta2_4 = list(map(float, delta2_4))
H12_4 = list(map(float, H12_4))

full_path5 = path5 + 'H12.dat'
data5 = open(full_path5)
x5, delta_5, delta1_5, delta2_5, H12_5 = arqui3.arquivo_dat(data5)
x5 = list(map(float, x5))
delta1_5 = list(map(float, delta1_5))
delta2_5 = list(map(float, delta2_5))
H12_5 = list(map(float, H12_5))

for i in range(len(x2)):
	x2[i] = x2[i]*100
for i in range(len(x3)):
	x3[i] = x3[i]*100
for i in range(len(x4)):
	x4[i] = x4[i]*100
for i in range(len(x5)):
	x5[i] = x5[i]*100


########################################################

#delta1 = np.zeros(len(x2))
#l = x2[len(x2)-1]

#for i in range(len(x2)):
#	delta1[i] = 0.34*(5.0*np.sqrt((nu*x2[i])/U))

########################################################

plt.figure()

plt.plot(x1, delta1, '.g', label='experiment')
plt.plot(x2, delta1_2,'b',label='blasius')
plt.plot(x3, delta1_3,'k',label='nssc')
plt.plot(x4, delta1_4,'r',label='grad_nssc')
plt.plot(x5, delta1_5,'m',label='grad003')

plt.plot(x1, delta2, '.g')
plt.plot(x2, delta2_2,'b')#,label='blasius')
plt.plot(x3, delta2_3,'k')#,label='nssc')
plt.plot(x4, delta2_4,'r')#,label='nscc')
plt.plot(x5, delta2_5,'m')#,label='nscc')

plt.plot(x1, H12, '.g')
plt.plot(x2, H12_2,'b')#,label='blasius')
plt.plot(x3, H12_3,'k')#,label='nssc')
plt.plot(x4, H12_4,'r')#,label='nscc')
plt.plot(x5, H12_5,'m')#,label='nscc')


plt.legend()
plt.xlim(0.0, 1200.0)
plt.ylim(0.0,3.0)
plt.show()
