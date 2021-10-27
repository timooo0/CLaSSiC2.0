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
maxEnergyIndex = np.argmax(energies > 10)
print(maxEnergyIndex)
fourierLength = int(param[fNum]["steps"] / 2)
sideLength = int(np.sqrt(param[fNum]["atoms"]))

if not os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat") and not os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat"):
    print("Calculating")
    start = time.time()
    latticePosition = helper.positionSquareLattice(sideLength, param[fNum])
    qScatter =  helper.qAll(sideLength)
    I_total, eIndex = helper.runTransform(latticePosition, x[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
    print(f'duration: {time.time()-start}')
else:
    print("Pulling fourier data from file")
    I_total = np.fromfile(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
    I_total = I_total.reshape((sideLength, sideLength, maxEnergyIndex))
    eIndex = np.fromfile(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat").astype(np.int32)


eIndex = eIndex[energies[eIndex].argsort()]
values, bins, _ = plt.hist(energies[eIndex], bins=144)
plt.xlabel("Energy [meV]")
plt.ylabel("Counts [-]")
plt.title("Histogram of the energies of the 2D spin waves")
eIndex = eIndex[np.where(values>=1)]

extent = [-np.pi, np.pi, -np.pi, np.pi]
fig, ax = plt.subplots(4, 4)
for i in range(4):
    for j in range(4):
        ax[i][j].imshow(I_total[:, :, int(eIndex[4 * i + 1*j])], extent=extent)
        ax[i][j].set_aspect(abs(extent[1] - extent[0]) / abs(extent[3] - extent[2]))
        ax[i][j].set_title(f'E = {energies[int(eIndex[4 * i + 1*j])]:.3f}')
        ax[i][j].get_xaxis().set_ticks([])
        ax[i][j].get_yaxis().set_ticks([])
# im = plt.imshow(I_total[:,:,eIndex], extent = extent)

# plt.plot(q[1:,0], helper.constants["J_to_meV"] * 4 * param[fNum]["J"]*3.5*(1-np.cos(q[1:,0])), 'w')


plt.show()
