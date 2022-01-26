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
            "geometry" : int(data[14]),
            "nUnitCells" : int(data[15]), 
            "anisotropyStrength" : data[16]
        }
        param.append(parameters)
        x.append(np.reshape(data[param[i]["offset"]+0::3], (param[i]["steps"], param[i]["atoms"])).T)
        y.append(np.reshape(data[param[i]["offset"]+1::3], (param[i]["steps"], param[i]["atoms"])).T)
        z.append(np.reshape(data[param[i]["offset"]+2::3], (param[i]["steps"], param[i]["atoms"])).T)
        
        # data for the different geometries:
        if parameters["geometry"] == 0:
            parameters["geometry"] = "single"
            parameters["nDimensions"] = 1
            parameters["nUnitCells"] = 1
            parameters["basisPosition"] = np.array([[0., 0., 0.]])
        elif parameters["geometry"] == 1:
            parameters["geometry"] = "line"
            parameters["nDimensions"] = 1
            parameters["basisPosition"] = np.array([[0., 0., 0.]])
            parameters["unitVectors"] = np.array([[1., 0., 0.]])
        elif parameters["geometry"] == 2:
            parameters["geometry"] = "square"
            parameters["nDimensions"] = 2
            parameters["basisPosition"] = np.array([[0., 0., 0.]])
            parameters["unitVectors"] = np.array([[1., 0., 0.], [0., 1., 0.]])
        elif parameters["geometry"] == 3:
            parameters["geometry"] = "triangle"
            parameters["nDimensions"] = 2
            parameters["basisPosition"] = np.array([[0., 0., 0.]])
            parameters["unitVectors"] = np.array([[1., 0., 0.], [np.cos(np.pi/3.), np.sin(np.pi/3.), 0.]])
        elif parameters["geometry"] == 4:
            parameters["geometry"] = "kagome"
            parameters["nDimensions"] = 2
            parameters["basisPosition"] = np.array([[0., 0., 0.], [1., 0., 0.], [1. * np.cos(np.pi/3.), 1. * np.sin(np.pi/3.), 0.]])
            parameters["unitVectors"] = np.array([[2., 0., 0.], [2. * np.cos(np.pi/3.), 2. * np.sin(np.pi/3.), 0.]])

        i += 1
        print(x[-1].shape)
        path = basePath[:-4] + str(i) + basePath[-4:]
    if i==0:
        print("Could not open the file!")
    
    x = np.stack(x, axis=0)
    y = np.stack(y, axis=0)
    z = np.stack(z, axis=0)


    return param, x, y, z

def scatterReal(q, position, spin, param):
    input = pyfftw.empty_aligned((param["atoms"], param["steps"]), dtype='float64')
    output = pyfftw.empty_aligned((param["atoms"], int(param["steps"]//2+1)), dtype='complex128')
    ic_sum = 0
    for i in range(param["atoms"]):

        input[i, :] = spin[i, :]
        ic_sum += np.exp(-1j * np.dot(q, position[i])) * spin[i, 0]

    output = pyfftw.FFTW(input, output, axes=(1,), flags=('FFTW_MEASURE',), threads=8)()
    for i in range(param["atoms"]):
        output[i,:] *= np.exp(1j * np.dot(q, position[i]))
    output = np.sum(output, axis=0)
    output *= ic_sum
    output = np.abs(output)
    return output

def scatter(q, position, spin, param):
    input = pyfftw.empty_aligned((param["atoms"], param["steps"]), dtype='complex128')
    output = pyfftw.empty_aligned((param["atoms"], param["steps"]), dtype='complex128')
    ic_sum = 0
    for i in range(param["atoms"]):

        input[i, :] = spin[i, :]
        ic_sum += np.exp(-1j * np.dot(q, position[i])) * spin[i, 0]

    output = pyfftw.FFTW(input, output, axes=(1,), flags=('FFTW_MEASURE',), threads=8)()
    for i in range(param["atoms"]):
        output[i,:] *= np.exp(1j * np.dot(q, position[i]))
    output = np.sum(output, axis=0)
    output *= ic_sum
    return output

def runTransform(latticePosition, spin, q, size, maxEnergyIndex, param):
    I_total = np.zeros((q.shape[0], maxEnergyIndex), dtype="complex128")
    for i in range(q.shape[0]):
        if i < q.shape[0]-1:
            print(f'progress: {i/q.shape[0]*100:.2f} %', end='\r')
        else:
            print(f'progress: finshed')
        I_aa = scatterReal(q[i, :], latticePosition, spin, param)
        # if i==1:
        #     frequencies = np.fft.fftfreq(param["steps"], param["dt"])
        #     energies = 4.1357e-12 * frequencies
        #     plt.plot(energies[:], np.abs(I_aa[:]), label="abs")
        #     plt.plot(energies[:], I_aa[:].real, label='real')
        #     plt.plot(energies[:], I_aa[:].imag, label='imag')
        I_total[i, :] = I_aa[:maxEnergyIndex]
    return I_total

def getPositions(fileNumber):
    positionBasePath = os.getcwd()+"\\CLaSSiC2.0\\data\\position.csv"
    positions = np.loadtxt(positionBasePath[:-4] + str(fileNumber) + positionBasePath[-4:], delimiter=", ")

    return positions

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

def scatterTriangle(size):
    q = np.array([0, 0, 0])
    print(f'triangluar lattice size: {size}')
    q = reciprocalPath([0,0,0], [np.pi, -1/np.sqrt(3)*np.pi, 0], q, size)
    q = reciprocalPath([np.pi,-1/np.sqrt(3)*np.pi,0], [np.pi, 0, 0], q, size)
    q = reciprocalPath([np.pi, 0, 0], [0,0,0], q, size)
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
        for i in range(4):
            ax.axvline((i)*length)
        ax.set_xticks([0, length, 2*length, 3*length])
        ax.set_xticklabels([r'$\Gamma$',r'$M$',r'$R$',r'$\Gamma$'])
    elif structure == "triangle":
        for i in range(4):
            ax.axvline((i)*length)
        ax.set_xticks([0, length, 2*length, 3*length])
        ax.set_xticklabels([r'$\Gamma$',r'$K$',r'$M$',r'$\Gamma$'])

def plotTheory(structure, length, param, ax, c):

    if structure=="line":
        q = scatterLine(int(length/2))
        xData = np.linspace(0,int((q.shape[0]))-1, q.shape[0])

        if param["J"] > 0:
            yData = np.abs(constants["J_to_meV"]*(4*param["J"]*3.5*(1-np.cos(q[:,0]))-constants["gFactor"]*constants["bohrMagneton"]*(param["magneticField"][-1]-param["anisotropyStrength"])))
        else:
            a1 = (4*param["J"]*3.5)**2*(1-np.cos(q[:,0])**2)
            a2 = -8*np.abs(param["J"])*3.5*constants["gFactor"]*constants["bohrMagneton"]*param["anisotropyStrength"]
            a3 = (constants["gFactor"]*constants["bohrMagneton"]*param["anisotropyStrength"])**2
            b = a1+a2+a3
            # yData = constants["J_to_meV"]*(np.sqrt((4*param["J"]*3.5)**2*(1-np.cos(q[:,0])**2)+(constants["gFactor"]*constants["bohrMagneton"]*param["anisotropyStrength"])**2)-constants["gFactor"]*constants["bohrMagneton"]*param["magneticField"][-1])
            yData = constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)-constants["gFactor"]*constants["bohrMagneton"]*param["magneticField"][-1])
            if param["magneticField"][-1] != 0:
                ax.plot(xData, yData, constants["colors"][c], label=f'B = {param["magneticField"][-1]:.0f} T, Anis = {param["anisotropyStrength"]} T', zorder=0)
                yData = constants["J_to_meV"]*(np.sqrt(a1 + a2 + a3)+constants["gFactor"]*constants["bohrMagneton"]*param["magneticField"][-1])

    if structure=="square":
        q = scatterSquare(int(length/2))
        xData = np.linspace(0,int((q.shape[0]))-1, q.shape[0])
        if param["J"] > 0:
             yData = np.abs(constants["J_to_meV"]*(8*param["J"]*3.5*(1-0.5*(np.cos(q[:,0])+np.cos(q[:,1])))-constants["gFactor"]*constants["bohrMagneton"]*(param["magneticField"][-1]-param["anisotropyStrength"])))

    if structure=="triangle":
        xData = []
        yData = []
    if param["magneticField"][-1] != 0:
            labelText = f'B = {param["magneticField"][-1]:.0f} T, Anis = {param["anisotropyStrength"]} T'
    else:
        labelText = f'Anis = {param["anisotropyStrength"]:.0f} T'

    ax.plot(xData, yData, constants["colors"][c], label=labelText)
    
