import numpy as np
import matplotlib.pyplot as plt


def arquivo_dat(data):
	
	x = []
	y = []
	vel_u = []
	vel_v = []
	vort_z = []
	line1 = data.readline()
	print line1
	line2 = data.readline()
	print line2
	sep = line2.split(' ')
	print sep
	n = len(sep)
	I = sep[10]
	print I
	J = sep[21]
	print J

	line = 'teste'
	while (line != ''):
		line = data.readline()
		if (line == ''):
			return I, J, x, y, vel_u, vel_v, vort_z
		
		line  = line.replace('D','E')
		sep = line.split(' ')
		n = len(sep)
		cont = 0
		for i in range(n):
			if (sep[i] != ''):
				if (cont == 0):
					x.append(sep[i])
					cont += 1

				elif (cont == 1):
					y.append(sep[i])
					cont += 1

				elif (cont == 2):
					vel_u.append(sep[i])
					cont += 1

				elif (cont == 3):
					vel_v.append(sep[i])
					cont += 1

				elif (cont == 4):
		 			var = sep[i]
					var = var.replace('\n','')
					vort_z.append(var)
					cont +=1

		cont = 0

	return I, J, x, y, vel_u, vel_v, vort_z	

def matrizar(I, J, M):
	
	A = np.zeros((J,I))
	cont = 0
	for j in range(J):
        	for i in range(I):
                	A[j,i] = M[cont]
                	cont +=1
	
	return A


path = '../SB/'
arq = 'basens.dat'
full_path = path + arq
print full_path
data = open(full_path)

I, J, x, y, vel_u, vel_v, vort_z = arquivo_dat(data)

print vort_z[0]
x = list(map(float, x))
y = list(map(float, y))
vel_u = list(map(np.float128, vel_u)) 
vel_v = list(map(np.float128, vel_v))
vort_z = list(map(np.float128, vort_z))

I = int(I)
J = int(J)

X = matrizar(I, J, x)
Y = matrizar(I, J, y)
U = matrizar(I, J, vel_u)
V = matrizar(I, J, vel_v)
Wz = matrizar(I, J, vort_z)

plt.figure()
plt.contourf(X, Y, U)
plt.xlabel('x')
plt.xlim(0,16)
plt.ylabel('y')
plt.title('U')
plt.show()
'''
plt.figure()
plt.contourf(X, Y, V)
plt.xlabel('x')
plt.ylabel('y')
plt.title('V')
plt.show()

plt.figure()
plt.contourf(X, Y, Wz)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Wz')
'''
plt.show()






