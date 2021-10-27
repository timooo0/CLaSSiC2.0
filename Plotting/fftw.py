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




def scatter(q, position):
    input = pyfftw.empty_aligned((param[fNum]["atoms"], param[fNum]["steps"]), dtype='complex128')
    output = pyfftw.empty_aligned((param[fNum]["atoms"], param[fNum]["steps"]), dtype='complex128')
    ic_sum = 0
    for i in range(param[fNum]["atoms"]):
        q_dot_lattice = np.exp(1j * np.dot(q, position[i]))

        input[i, :] = q_dot_lattice * x[fNum, i, :]

        ic_sum += cmath.exp(-1j * np.dot(q, position[i])) * x[fNum, i, 0]


    output = pyfftw.FFTW(input, output, axes=(1,), flags=('FFTW_MEASURE',), threads=8)()
    output = np.sum(output, axis=0)
    output *= ic_sum
    output = np.power(np.abs(output), 2)
    return output[:fourierLength].real


def runTransform():
    lattice_position = []
    for i in range(sideLength):
        for j in range(sideLength):
            lattice_position.append([j, i, 0])
    lattice_position = np.array(lattice_position).reshape(param[fNum]["atoms"], 3)

    I_total = np.zeros((sideLength, sideLength, maxEnergyIndex))

    q_x = np.array([0, 0, 0])
    eIndex = np.array([], dtype=np.float32)
    for i in range(sideLength):
        print(f"q_x: {i}")
        for j in range(sideLength):
            q_x = np.vstack((q_x, np.array([np.pi * (2 * i / sideLength - 1), np.pi * (2 * j / sideLength - 1), 0])))
            I_aa = scatter(q_x[-1, :], lattice_position)
            eIndex = np.append(eIndex, np.argmax(I_aa))
            I_total[i, j, :] = I_aa[:maxEnergyIndex]

    I_total.tofile(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
    eIndex.tofile(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat")

    return I_total, eIndex.astype(np.int32)


if not os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat") and not os.path.exists(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat"):
    print("Recalculating")
    start = time.time()
    I_total, eIndex = runTransform()
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
