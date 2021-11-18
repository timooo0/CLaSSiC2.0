import os
import numpy as np
import pyfftw
import cmath
import matplotlib.pyplot as plt


constants = {
    "hbar" : 1.054571817e-34,
    "boltzmann" : 1.38064852e-23,
	"gFactor" : -2.002,
	"bohrMagneton" : 9.274009994e-24,
	"gamma" : 1.760859644e11,
	
    "J_to_meV" : 6.2415093433e21,
    "Hz_to_meV": 4.135667696e-12,

    "colors" : ['b','g','r','c','m','y','b'],

    "pathFourier" : os.getcwd() + "\\CLaSSiC2.0\\data\\fourier.dat",
    "pathEnergyIndex" : os.getcwd() + "\\CLaSSiC2.0\\data\\energyIndex.dat"
}

def getPath(path, i):
    return path[:-4] + str(i) + path[-4:]

def getData():
    param = []
    x, y, z = [], [], []
    basePath = os.getcwd()+"\\CLaSSiC2.0\\data\\data.dat"
    i = 0
    path = basePath[:-4] + str(i) + basePath[-4:]
    print(path)
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
            "length" : int(data[13]), 
            "anisotropyStrength" : data[14]
        }
        param.append(parameters)
        x.append(np.reshape(data[param[i]["offset"]+0::3], (param[i]["steps"], param[i]["atoms"])).T)
        y.append(np.reshape(data[param[i]["offset"]+1::3], (param[i]["steps"], param[i]["atoms"])).T)
        z.append(np.reshape(data[param[i]["offset"]+2::3], (param[i]["steps"], param[i]["atoms"])).T)

        i += 1
        print(x[-1].shape)
        path = basePath[:-4] + str(i) + basePath[-4:]
    if i==0:
        print("Could not open the file!")
    
    x = np.stack(x, axis=0)
    y = np.stack(y, axis=0)
    z = np.stack(z, axis=0)

    return param, x, y, z

def scatter(q, position, spin, param):
    input = pyfftw.empty_aligned((param["atoms"], param["steps"]), dtype='float64')
    output = pyfftw.empty_aligned((param["atoms"], int(param["steps"]//2+1)), dtype='complex128')
    ic_sum = 0
    for i in range(param["atoms"]):
        input[i, :] = spin[i, :]
        ic_sum += np.exp(-1j * np.dot(q, position[i])) * spin[i, 0]

    output = pyfftw.FFTW(input, output, axes=(1,), flags=('FFTW_MEASURE',), threads=8)()
    for i in range(param["atoms"]):
        # if i>30: print(np.sum(np.exp(1j * np.dot(q, position[i]))))
        output[i,:] *= np.exp(1j * np.dot(q, position[i]))
        # plt.plot(output[i,:])
    output = np.sum(output, axis=0)
    # for i in range(param["atoms"]):
    #     output[i] *= np.exp(-1j * np.dot(q, position[i])) * spin[i, 0]
    output *= ic_sum
    output = np.abs(output)
    return output

def runTransform(latticePosition, spin, q, size, maxEnergyIndex, param):
    I_total = np.zeros((q.shape[0], maxEnergyIndex))
    for i in range(q.shape[0]):
        print(i)
        I_aa = scatter(q[i, :], latticePosition, spin, param)
        I_total[i, :] = I_aa[:maxEnergyIndex]

    # for i in range(q.shape[0]):
    #     I_total[i, :] = I_total[i, :]/np.max(I_total[i, :])
    print(I_aa.shape, maxEnergyIndex)
    return I_total

def positionLine(size, param):
    latticePosition = []
    for i in range(size):
        latticePosition.append([i, 0, 0])
    latticePosition = np.array(latticePosition).reshape(param["atoms"], 3)

    return latticePosition

def positionSquare(size, param):
    latticePosition = []
    for i in range(size):
        for j in range(size):
            latticePosition.append([j, i, 0])
    latticePosition = np.array(latticePosition).reshape(param["atoms"], 3)

    return latticePosition

def scatterLine(size):
    q = np.array([0, 0, 0])
    q = reciprocalPath([-np.pi, 0, 0], [0, 0, 0], q, size)
    q = reciprocalPath([0, 0, 0], [np.pi, 0, 0], q, size)
    q = np.delete(q, 0, axis=0)
    return q

def scatterSquare(size):
    q = np.array([0, 0, 0])
    q = reciprocalPath([0,0,0], [np.pi, 0, 0], q, size)
    q = reciprocalPath([np.pi,0,0], [np.pi, np.pi, 0], q, size)
    q = reciprocalPath([np.pi,np.pi,0], [0, 0, 0], q, size)
    q = np.delete(q, 0, axis=0)

    return q

def scatterAll(size):
    q = np.array([0, 0, 0])
    for i in range(size):
        for j in range(size):
            q = np.vstack((q, np.array([np.pi * (2 * i / size - 1), np.pi * (2 * j / size - 1), 0])))

    q = np.delete(q, 0, axis=0)

    return q

def reciprocalPath(start, end, q, size):
    length = np.array(end) - np.array(start)
    for i in range(size): 
        q = np.vstack((q, np.array([start[0] + length[0]* i / size, start[1] + length[1]* i / size, start[2] + length[2]* i / size])))
    return q

def plotMarker(structure, length, ax):
    length = int(length/2)
    if structure == "line":
        ax.axvline(length)
        ax.set_xticks([0, length, 2*length])
        ax.set_xticklabels([r'$R$',r'$\Gamma$',r'$R$'])
    elif structure == "square":
        for i in range(3-1):
            ax.axvline((1+i)*length)
        ax.set_xticks([0, length, 2*length, 3*length])
        ax.set_xticklabels([r'$\Gamma$',r'$M$',r'$R$',r'$\Gamma$'])

def plotTheory(structure, length, param, ax, c):

    if structure=="line":
        q = scatterLine(int(length/2))
        xData = np.linspace(0,int((q.shape[0]))-1, q.shape[0])
        if param["J"] > 0:
            yData = np.abs(constants["J_to_meV"]*(4*param["J"]*3.5*(1-np.cos(q[:,0]))-constants["gFactor"]*constants["bohrMagneton"]*(param["magneticField"][-1]-param["anisotropyStrength"])))
        else:
            a1 = (4*param["J"]*3.5)**2*(1-np.cos(q[:,0])**2)
            a2 = 8*param["J"]*3.5**2*constants["gFactor"]*constants["bohrMagneton"]*param["anisotropyStrength"]/(2)
            a3 = 3.5**2*(constants["gFactor"]*constants["bohrMagneton"]*param["anisotropyStrength"]/2)**2
            b = a1+a3
            # yData = constants["J_to_meV"]*(np.sqrt((4*param["J"]*3.5)**2*(1-np.cos(q[:,0])**2)+(constants["gFactor"]*constants["bohrMagneton"]*param["anisotropyStrength"])**2)-constants["gFactor"]*constants["bohrMagneton"]*param["magneticField"][-1])
            yData = constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)-constants["gFactor"]*constants["bohrMagneton"]*param["magneticField"][-1])
            ax.plot(xData, yData, constants["colors"][c], label=f'B = {param["magneticField"][-1]:.0f} T, Anis = {param["anisotropyStrength"]} T', zorder=0)
            yData = constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)+constants["gFactor"]*constants["bohrMagneton"]*param["magneticField"][-1])
        ax.plot(xData, yData, constants["colors"][c], label=f'B = {param["magneticField"][-1]:.0f} T, Anis = {param["anisotropyStrength"]} T', zorder=0)