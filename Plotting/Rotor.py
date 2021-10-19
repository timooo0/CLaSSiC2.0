import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()

energy, angle = [], []
print(x.shape)

for fNum in range(x.shape[0]):
    energy.append(helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"])[np.argmax(np.fft.fft(x[fNum, 0, :]))])
    angle.append(90-np.arccos(z[fNum, 0, -1])*180/np.pi)

theoryAngle  = np.linspace(0,90,100)
theoryEnergy = helper.constants["J_to_meV"] * 4 * param[fNum]["J"] * 7/2 * np.sin(np.radians(theoryAngle))


fig, ax = plt.subplots()
ax.plot(theoryAngle, theoryEnergy, 'r')
ax.plot(angle, np.abs(energy),'bo')
ax.set_xlabel(r'$\theta$ [degrees]')
ax.set_ylabel(r'Energy [meV]')
ax.set_title('Rotor')

# time = np.linspace(0,param[fNum][fNum]["dt"]*N,N)
# plt.plot(time,x[:, 0], 'r')
# ax[1].plot(4.135665538536e-12*np.real(np.fft.fftfreq(np.size(x),d=dt)), np.fft.fft(x))
# ax[1].set_xlim(-1,1)
#
# # plt.axis('equal')
plt.show()