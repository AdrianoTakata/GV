import numpy as np
import matplotlib.pyplot as plt

N = 100
R_lambda = 1145.0
lambda_z = 0.023

x_hat = np.linspace(0.001, 0.4, N)

x = x_hat * R_lambda * (lambda_z/2*np.pi)
print x

plt.figure()
plt.plot(x)
plt.show()
