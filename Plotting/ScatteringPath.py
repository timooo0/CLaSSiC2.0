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
print(f'energies: {energies}, size: {energies.shape}, max: {np.max(energies)}')
maxEnergyIndex = np.argmax(energies > 10)
maxEnergyIndex = int(param[fNum]["steps"]//2+1)
print(param[fNum])
fourierLength = int(param[fNum]["steps"] / 2)

print(f'z sum: {np.mean(np.sum(z[fNum],axis=0))}')
if os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier0.dat"):
    os.remove(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier0.dat")

def getTransform(structure):
    if structure == "line":
        sideLength = int(param[fNum]["atoms"])
        latticePosition = helper.positionLine(sideLength, param[fNum])
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/2)
        qScatter = helper.scatterLine(scatterLength)
    elif structure == "square":
        sideLength = int(np.sqrt(param[fNum]["atoms"]))
        latticePosition = helper.positionSquare(sideLength, param[fNum])
        scatterLength = int(sideLength/2) if param[fNum]["J"] >= 0 else int(sideLength/4)
        qScatter = helper.scatterSquare(scatterLength)

    nScatter = qScatter.shape[0]
    pathFourier = helper.getPath(helper.constants["pathFourier"], fNum)
    pathEnergyIndex = helper.getPath(helper.constants["pathEnergyIndex"], fNum)
    if not os.path.exists(pathFourier) and not os.path.exists(pathEnergyIndex):
        print("Calculating")
        start = time.time()
        I_total_x = helper.runTransform(latticePosition, x[fNum,:,:], qScatter, sideLength, maxEnergyIndex, param[fNum])
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
point = 0
while fNum < x.shape[0]:
    I_total, nScatter, sideLength = getTransform("line")
    xData = np.linspace(0,nScatter-1,nScatter)
    yData, yData2, zData = [], [], []
    maxF = np.log10(np.max(I_total))
    for j in range(nScatter):
        yData.append(energies[np.argmax(I_total[j,:])])
        zData.append(max([np.log10(np.max(I_total[j,:]))-maxF, -5]))
        # ax[0].plot(xData[j],energies[np.argmax(I_total[j,:])], c=plt.get_cmap("autumn")((-np.log10(np.max(I_total[j,:]))-maxF)/3),marker='o')
        # print(np.max(I_total[j,:])/maxF, (np.log10(np.max(I_total[j,:]))-maxF)/3)
        # yData2.append(energies[np.argsort(I_total[j,:])[-2]])
        # peaks, _ = find_peaks(I_total[j, :], height=0.5, distance=50)
        # print(f'j: {j}, size: {peaks.shape}')
        # if peaks.shape[0] == 0:
        #     yData.append(energies[np.argmax(I_total[j,:])])
        #     yData2.append(energies[np.argmax(I_total[j,:])])
        # elif peaks.shape[0] == 1:
        #     yData.append(energies[peaks[0]])
        #     yData2.append(energies[peaks[0]])
        # else:
        #     yData.append(energies[peaks[0]])
        #     yData2.append(energies[peaks[1]])

    # ax[0].plot(xData,yData2, helper.constants["colors"][fNum+1]+'o')
    # ax[0].plot(xData,yData3, helper.constants["colors"][fNum]+'o')
    helper.plotTheory('line',sideLength, param[fNum], ax[0], fNum)
    im = ax[0].scatter(xData,yData,c=zData,cmap="Wistia", zorder=10)
    fig.colorbar(im, ax=ax[0])

    fNum = fNum + 1

# ax[1].plot(energies, I_total[point,:].real)
ax[1].plot(energies[:maxEnergyIndex-1], I_total[point,:-1])
# ax[1].set_yscale("log")
peaks, _ = find_peaks(I_total[point, :], height=0.5, distance=250)
# ax[1].plot([energies[peaks[0]], energies[peaks[1]]], [I_total[point, peaks[0]], I_total[point, peaks[1]]], 'x')
ax[0].axvline(xData[point])
fNum = fNum - 1
ax[0].set_ylabel('Energy [meV]')
ax[0].legend(loc="lower right")
if param[fNum]["J"] > 0:
    ax[0].set_title('Ferromagnetic spin chain')
else:
    ax[0].set_title('Antiferromagnetic spin chain')
helper.plotMarker('line', sideLength, ax[0])

plt.show()
