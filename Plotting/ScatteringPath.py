from numpy.core.fromnumeric import argmax
import helper
import numpy as np
import os
import time
import cmath
import pyfftw
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = np.argmax(energies > 15)
print(maxEnergyIndex)
fourierLength = int(param[fNum]["steps"] / 2)
sideLength = int(np.sqrt(param[fNum]["atoms"]))

os.remove(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
os.remove(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat")
if not os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat") and not os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat"):
    print("Calculating")
    start = time.time()
    latticePosition = helper.positionSquareLattice(sideLength, param[fNum])
    qScatter =  helper.qSquare(sideLength)
    I_total, eIndex = helper.runTransform(latticePosition, y[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
    I_total = I_total.reshape((qScatter.shape[0], maxEnergyIndex))
    print(f'duration: {time.time()-start}')
else:
    print("Pulling fourier data from file")
    I_total = np.fromfile(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
    I_total = I_total.reshape((3*sideLength, maxEnergyIndex))
    eIndex = np.fromfile(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat").astype(np.int32)


print(I_total.shape)
# extent = [0, 36, 0, 6]
# plt.imshow(I_total, extent=extent)
xData = np.linspace(0,3*sideLength,3*sideLength)
yData = []
for i in range(3*sideLength):
    yData.append(energies[argmax(I_total[i,:])])
plt.plot(xData,yData, 'bo')
for i in range(3-1):
    plt.axvline((1+i)*sideLength)
plt.xticks([0, sideLength, 2*sideLength, 3*sideLength], [r'$\Gamma$',r'$M$',r'$R$',r'$\Gamma$'])

plt.show()
