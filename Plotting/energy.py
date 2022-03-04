import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0

fig, ax = plt.subplots()
end = int(1e4)
time = np.linspace(0,param[fNum]["dt"]*param[fNum]["steps"],param[fNum]["steps"])
energy = helper.getEnergy(fNum)
ax.plot(time[:end], energy[:end])
# ax.set_yscale("log")


plt.show()
