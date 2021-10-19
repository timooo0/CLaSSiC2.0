import helper
import numpy as np
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum  = 0
atom = 0

print(x.shape)
fig, ax = plt.subplots(1,3)
time = np.linspace(0,param[fNum]["steps"]*param[fNum]["dt"],param[fNum]["steps"])
#energy = data[3::N]


print(f'z difference: {np.max(z[fNum, atom, :])-np.min(z[fNum, atom, :])}')
print(f'final length: {np.sqrt(x[fNum, atom, -1]**2 + y[fNum, atom, -1]**2 + z[fNum, atom, -1]**2)}')
precession = helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"])[np.argmax(np.fft.fft(x[fNum, atom, :]))]
theory = helper.constants["Hz_to_meV"]*helper.constants["gamma"]*param[fNum]["magneticField"][-1]/(2*np.pi)
print(f'precession energy: {precession} meV, theoretical: {theory} meV, error: {theory/precession}')

end = int(1e5)
ax[0].plot(time[:end], x[fNum, atom, :end],'r')
ax[0].plot(time[:end], y[fNum, atom, :end],'b')
# ax[0].plot(time[:end], z[fNum, atom, :end],'g')
ax[0].set_xlabel('time [s]')
ax[0].set_ylabel('spin [-]')
ax[0].set_title('Spin evolution')

ax[1].plot(helper.constants["Hz_to_meV"]*np.fft.fftfreq(param[fNum]["steps"],d=param[fNum]["dt"]), np.fft.fft(x[fNum, atom, :]).real)
ax[1].set_xlim(-2, 2)
ax[1].set_xlabel('Energy [meV]')
ax[1].set_ylabel('relative count [-]')
ax[1].set_title('Fourier transform')


ax[2].plot(x[fNum, atom, :end], y[fNum, atom, :end],'r')
ax[2].set_xlabel('x')
ax[2].set_xlabel('y')
ax[2].set_title('x-y plane')
# # plt.axis('equal')
plt.show()
