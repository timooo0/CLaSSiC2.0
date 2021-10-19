import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0
print(x[fNum, :,:].shape)
# print(data[offset:offset+36])
# for i in range(12):
#     print(f'x: {x[0,i]}, y: {y[0, i]}')
# print(f'z difference: {z[:, -1]-z[:, -1]}')
fig, ax = plt.subplots(1,2)
end = int(1e4)
time = np.linspace(0,param[fNum]["dt"]*param[fNum]["steps"],param[fNum]["steps"])
for i in range(1):
    for j in range(2):
        print(3*i+j, x[0, 3*i+j])
        ax[j].plot(time[:end],x[fNum, 3*i+j, :end], 'r')
        ax[j].plot(time[:end],y[fNum, 3*i+j, :end], 'b')
        ax[j].plot(time[:end],z[fNum, 3*i+j, :end], 'g')


plt.show()
