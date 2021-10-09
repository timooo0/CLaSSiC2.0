import os
import numpy as np
import matplotlib.pyplot as plt

data = np.fromfile(os.getcwd()+"\\CLaSSiC2.0\\data\\data0.dat", dtype=np.double)

offset = int(data[0])
dt = data[1]
J = data[2]
Lambda = data[3]
magneticField = data[4:7]
anisotropy = data[7:10]
temperature = data[10]
length = int(data[11])
print(f'offset: {offset}')
print(f'dt: {dt}')
print(f'J: {J}')
print(f'lambda: {Lambda}')
print(f'magnetic field: {magneticField}')
print(f'anisotropy: {anisotropy}')
print(f'temperature: {temperature}')
print(f'length: {length}')

x = data[offset::length]
y = data[offset+1::length]
z = data[offset+2::length]




fig, ax = plt.subplots(1,3)
#


time = np.linspace(0,np.size(x)*dt,np.size(x))
#energy = data[3::N]

print(f'size of x: {np.size(x)}')
print(f'size of spin: {np.size(data)-offset}')
print(f'z difference: {np.max(z)-np.min(z)}')
print(np.size(np.fft.fft(x)))
print('plotting...')
print(f'precession energy: {4.135665538536e-12*np.fft.fftfreq(np.size(x),d=dt)[np.argmax(np.fft.fft(x))]} meV')

end = 1000_00
ax[0].plot(time[:], x[:],'r')
ax[0].plot(time[:], y[:],'b')
ax[0].plot(time[:], z[:],'g')
ax[0].set_xlabel('time [s]')
ax[0].set_ylabel('spin [-]')
ax[0].set_title('Spin evolution')

ax[1].plot(4.135665538536e-12*np.fft.fftfreq(np.size(x),d=dt), np.fft.fft(x).real/np.size(x))
ax[1].set_xlim(-2, 2)
ax[1].set_xlabel('Energy [meV]')
ax[1].set_ylabel('relative count [-]')
ax[1].set_title('Fourier transform')


ax[2].plot(x[:end], y[:end],'r')
ax[2].set_xlabel('x')
ax[2].set_xlabel('y')
ax[2].set_title('x-y plane')
# # plt.axis('equal')
plt.show()
