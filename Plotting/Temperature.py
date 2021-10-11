import os
import numpy as np
import matplotlib.pyplot as plt
basePath = os.getcwd()+"\\CLaSSiC2.0\\data\\data.dat"
i = 0
path = basePath[:-4] + str(i) + basePath[-4:]


succes = True
temperatureArray = []
spinZ = []

while(succes):
    try:
        with open(path, "rb") as f:
            data = np.fromfile(f, dtype=np.double)

        offset = int(data[0])
        temperatureArray.append(data[10])
        length = int(data[11])
        magneticField = data[4:7]
        z = data[offset + 2::length]
        spinZ.append(np.mean(z))

        i+=1
        path = basePath[:-4] + str(i) + basePath[-4:]
        print(f'file {i}')
    except IOError:
        if i==0:
            print('Error While Opening the file!')
        else:
            succes = False

#


print('plotting...')
fig, ax = plt.subplots()

y = 2.002 * 7/2 * magneticField[-1]*9.274009994e-24/1.38064852e-23

temperature = np.linspace(1,100,50)
y = y/np.array(temperature)
langevin = np.cosh(y)/np.sinh(y) - 1/y
print(f'Magnetic field: {magneticField}')

ax.plot(temperature, langevin,'b')
ax.plot(temperatureArray, np.abs(spinZ), 'r.')
ax.set_xlabel('$T [K]$')
ax.set_ylabel('$S [-]$')
#
# # plt.axis('equal')
plt.show()
f.close()