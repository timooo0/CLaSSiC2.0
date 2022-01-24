import helper
import numpy as np
from scipy.signal import find_peaks
import os
import time
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
maxEnergyIndex = int(param[fNum]["steps"]//2+1)
fourierLength = int(param[fNum]["steps"] / 2)

print(f'z sum: {np.mean(np.sum(z[fNum],axis=0))}')
if os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier0.dat"):
    os.remove(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier0.dat")

def getTransform(structure):
    if structure == "line":
        sideLength = int(param[fNum]["atoms"])
        latticePosition = helper.positionLine(sideLength)
        scatterLength = int(sideLength/2)
        qScatter = helper.scatterLine(scatterLength)

    elif structure == "line2":
        sideLength = int(param[fNum]["atoms"]/2)
        latticePosition = helper.positionLine(sideLength)
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/2)
        qScatter = helper.scatterLine(scatterLength)

    elif structure == "square":
        sideLength = int(np.sqrt(param[fNum]["atoms"]))
        latticePosition = helper.positionSquare(sideLength)
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/4)
        qScatter = helper.scatterSquare(scatterLength)

    elif structure == "triangle":
        sideLength = int(np.sqrt(param[fNum]["atoms"]))
        latticePosition = helper.positionTriangle(sideLength)
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/2)
        qScatter = helper.scatterTriangle(scatterLength)

    nScatter = qScatter.shape[0]
    print(f'nScatter: {qScatter.shape[0]}')
    pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
    pathEnergyIndex = helper.getPath(helper.constants["pathEnergyIndex"], fNum)
    if not os.path.exists(pathFourier) and not os.path.exists(pathEnergyIndex):
        print("Calculating")
        start = time.time()
        print("x: ")
        I_total_x = helper.runTransform(latticePosition, x[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
        print("y: ")
        I_total_y = helper.runTransform(latticePosition, y[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
        # I_total_z = helper.runTransform(latticePosition, z[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
        I_total = I_total_x**2 + I_total_y**2 #+ I_total_z**2
        I_total.tofile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
        print(f'duration: {time.time()-start}')
    else:
        print("Pulling fourier data from file")
        I_total = np.fromfile(pathFourier)
        I_total = I_total.reshape((nScatter, maxEnergyIndex))
    
    return I_total, nScatter, sideLength

fig, ax = plt.subplots(2)
point = 13
structure = "triangle"
while fNum < x.shape[0]:
    I_total, nScatter, sideLength = getTransform(structure)
    xData = np.linspace(0,nScatter-1,nScatter)
    yData, zData = [], []
    maxF = np.log10(np.max(I_total))
    for j in range(nScatter):
        yData.append(np.abs(energies[np.argmax(I_total[j,:])]))
        zData.append(max([np.log10(np.max(I_total[j,:]))-maxF, -5]))

    # ax[0].plot(xData,yData, helper.constants["colors"][fNum+1]+'o')
    # ax[0].plot(xData,yData3, helper.constants["colors"][fNum]+'o')
    print(xData.dtype, type(yData[0]), type(zData))
    helper.plotTheory(structure, sideLength, param[fNum], ax[0], fNum)
    im = ax[0].scatter(xData,yData,c=zData,cmap="Wistia", zorder=10)
    fig.colorbar(im, ax=ax[0])

    fNum = fNum + 1

fNum = fNum - 1

ax[1].plot(energies[:maxEnergyIndex-1], I_total[point,:maxEnergyIndex-1])
# ax[1].set_yscale("log")
ax[0].axvline(xData[point], c='r')
ax[0].set_ylabel('Energy [meV]')
ax[0].legend(loc="lower right")

if param[fNum]["J"] > 0:
    ax[0].set_title('Ferromagnetic spin chain')
else:
    ax[0].set_title('Antiferromagnetic spin chain')
helper.plotMarker(structure, sideLength, ax[0])

plt.show()
