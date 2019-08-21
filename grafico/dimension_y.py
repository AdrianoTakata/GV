import numpy as np


dy = 5.0e-4
N = 297
y = 0.0
stf = 1.01

for i in range(1,N):
	y = y + dy*stf**i

print y
