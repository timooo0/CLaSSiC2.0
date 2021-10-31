import os
import numpy as np
import matplotlib.pyplot as plt
import helper

param, x, y, z = helper.getData()
temperatureArray = []
zArray = []
for i in range(x.shape[0]):
    temperatureArray.append(param[i]["temperature"])
    zArray.append(np.mean(z[i,0,int(param[i]["steps"]/2):]))

#


print('plotting...')
fig, ax = plt.subplots()

y = 2.002 * 7/2 * param[0]["magneticField"][-1]*9.274009994e-24/1.38064852e-23

temperature = np.linspace(1,param[-1]["temperature"],50)
y = y/np.array(temperature)
langevin = np.cosh(y)/np.sinh(y) - 1/y
print(f'Magnetic field: {param[0]["magneticField"][-1]}')

ax.plot(temperature, langevin,'b')
ax.plot(temperatureArray, zArray, 'r.')
ax.set_xlabel('$T [K]$')
ax.set_ylabel('$S [-]$')
#
# # plt.axis('equal')
plt.show()