import os
import numpy as np

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
            "dt" : data[2],
            "steps" : int(data[3]),
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