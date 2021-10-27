import os
import numpy as np
import pyfftw
import cmath


constants = {
    "hbar" : 1.054571817e-34,
    "boltzmann" : 1.38064852e-23,
	"gFactor" : -2.002,
	"bohrMagneton" : 9.274009994e-24,
	"gamma" : 1.760859644e11,
	
    "J_to_meV" : 6.2415093433e21,
    "Hz_to_meV": 4.135667696e-12,
}

def getData():
    param = []
    x, y, z = [], [], []
    basePath = os.getcwd()+"\\CLaSSiC2.0\\data\\data.dat"
    i = 0
    path = basePath[:-4] + str(i) + basePath[-4:]
    while (os.path.exists(path)):
        data = np.fromfile(path, dtype=np.double)

        parameters = {
            "offset" : int(data[0]),
            "atoms" : int(data[1]),
            "dt" : data[2]*100,
            "steps" : int(data[3]/100),
            "J" : data[4],
            "lambda" : data[5],
            "magneticField" : data[6:9],
            "anisotropy" : data[9:12],
            "temperature" : data[12],
            "length" : int(data[13]) 
        }
        param.append(parameters)

        x.append(np.reshape(data[param[i]["offset"]+0::3], (param[i]["steps"], param[i]["atoms"])).T)
        y.append(np.reshape(data[param[i]["offset"]+1::3], (param[i]["steps"], param[i]["atoms"])).T)
        z.append(np.reshape(data[param[i]["offset"]+2::3], (param[i]["steps"], param[i]["atoms"])).T)

        i += 1
        path = basePath[:-4] + str(i) + basePath[-4:]
    
    x = np.stack(x, axis=0)
    y = np.stack(y, axis=0)
    z = np.stack(z, axis=0)

    return param, x, y, z

def scatter(q, position, spin, param):
    input = pyfftw.empty_aligned((param["atoms"], param["steps"]), dtype='complex128')
    output = pyfftw.empty_aligned((param["atoms"], param["steps"]), dtype='complex128')
    ic_sum = 0
    for i in range(param["atoms"]):
        q_dot_lattice = np.exp(1j * np.dot(q, position[i]))

        input[i, :] = q_dot_lattice * spin[i, :]

        ic_sum += cmath.exp(-1j * np.dot(q, position[i])) * spin[i, 0]


    output = pyfftw.FFTW(input, output, axes=(1,), flags=('FFTW_MEASURE',), threads=8)()
    output = np.sum(output, axis=0)
    output *= ic_sum
    output = np.power(np.abs(output), 2)
    return output[:int(param["steps"]/2)].real

def runTransform(latticePosition, spin, q, size, maxEnergyIndex, param):
    I_total = np.zeros((size, size, maxEnergyIndex))

    
    eIndex = np.array([], dtype=np.float32)
    for i in range(size):
        print(f"q_x: {i}")
        for j in range(size):
            I_aa = scatter(q[size*i+j, :], latticePosition, spin, param)
            eIndex = np.append(eIndex, np.argmax(I_aa))
            I_total[i, j, :] = I_aa[:maxEnergyIndex]

    I_total.tofile(os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat")
    eIndex.tofile(os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat")

    return I_total, eIndex.astype(np.int32)

def positionSquareLattice(size, param):
    latticePosition = []
    for i in range(size):
        for j in range(size):
            latticePosition.append([j, i, 0])
    latticePosition = np.array(latticePosition).reshape(param["atoms"], 3)

    return latticePosition

def qAll(size):
    q = np.array([0, 0, 0])
    for i in range(size):
        for j in range(size):
            q = np.vstack((q, np.array([np.pi * (2 * i / size - 1), np.pi * (2 * j / size - 1), 0])))

    q = np.delete(q, 0, axis=0)

    return q