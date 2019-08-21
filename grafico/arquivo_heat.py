import numpy as np

def arquivo_dat(data):

        var1 = []
        var2 = []
        var3 = []
        var4 = []
	var5 = []
	var6 = []
	var7 = []
	var8 = []
        line1 = data.readline()
        line2 = data.readline()
#	print 'start here'
        line = 'teste'
        while (line != ''):
                line = data.readline()
                if (line == ''):
			
			var1 = list(map(float, var1))
			var2 = list(map(float, var2)) 
			var3 = list(map(float, var3))
			var4 = list(map(float, var4))
			var5 = list(map(float, var5))
			var6 = list(map(float, var6))
			var7 = list(map(float, var7))
			var8 = list(map(float, var8))

                        return  var1, var2, var3, var4, var5, var6, var7, var8

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

				elif (cont == 4):
					var5.append(sep[i])
					cont += 1
				
				elif (cont == 5):
					var6.append(sep[i])
					cont += 1

				elif (cont == 6):
					var7.append(sep[i])
					cont += 1

				elif (cont == 7):
					var8.append(sep[i])
					cont += 1

                cont = 0

	var1 = list(map(float, var1))
	var2 = list(map(float, var2)) 
	var3 = list(map(float, var3))
	var4 = list(map(float, var4))
	var5 = list(map(float, var5))
	var6 = list(map(float, var6))
	var7 = list(map(float, var7))
	var8 = list(map(float, var8))

        return var1, var2, var3, var4, var5, var6, var7, var8

