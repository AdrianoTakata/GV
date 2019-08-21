import numpy as np

def arquivo_dat(data):

        var1 = []
        var2 = []
        var3 = []
	var4 = []
        line1 = data.readline()
        print line1
        line2 = data.readline()
        print line2
        sep = line2.split(' ')
	print sep
	I = sep[10]
        print 'i=', I
        J = sep[21]
        print 'j=', J
	K = sep[33]
	print 'k=', K

        line = 'teste'
        while (line != ''):
                line = data.readline()
                if (line == ''):
                        return  I, J, K, var1, var2, var3, var4

                line  = line.replace('D','E')
                sep = line.split(' ')
                n = len(sep)
                cont = 0
                for i in range(n):
                        if (sep[i] != ''):
                                if (cont == 0):
                                        var1.append(sep[i])
                                        cont += 1

                                elif (cont == 1):
                                        var2.append(sep[i])
                                        cont += 1

                                elif (cont == 2):
                                        var3.append(sep[i])
                                        cont += 1
				elif (cont == 3):
                                        var4.append(sep[i])
                                        cont += 1

                cont = 0

        return I, J, K, var1, var2, var3, var4


