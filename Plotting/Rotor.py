import numpy as np
import matplotlib.pyplot as plt
basePath = "C:/Users/timov/source/repos/CLaSSiC2.0/CLaSSiC2.0/data.dat"
i = 0
path = basePath[:-4] + str(i) + basePath[-4:]

success = True
angleArray = []
energyArray = []
print(path)
while(success):
    try:
        with open(path, "rb") as f:
            data = np.fromfile(f, dtype=np.double)

        offset = int(data[0])
        dt = data[1]
        J = data[2]
        length = int(data[11])
        x = data[offset::length]
        z = data[offset + 2::length]
        angleArray.append(np.degrees(np.arccos(z[0])))
        energyArray.append(6.241509e21*1.0545718e-34*2*np.pi*np.fft.fftfreq(np.size(x), d=dt)[np.argmax(np.fft.fft(x))])
        i+=1
        path = basePath[:-4] + str(i) + basePath[-4:]
        print(f'file {i}')
    except IOError:
        if i==0:
           print('Error While Opening the file!')
        else:
           success = False


fig, ax = plt.subplots()

print(f'size: {np.size(x)}')
print(f'z: {angleArray}')
print('plotting...')
end = 5000

theoryAngle  = np.linspace(0,90,100)
theoryEnergy = 6.241509e21 * 4 * J * 7/2 * np.cos(np.radians(theoryAngle))
ax.plot(theoryAngle, theoryEnergy, 'r')
ax.plot(angleArray[:end], np.abs(energyArray[:end]),'bo')
ax.set_xlabel(r'$\theta$ [degrees]')
ax.set_ylabel(r'Energy [meV]')
ax.set_title('Rotor')

# ax[1].plot(4.135665538536e-12*np.real(np.fft.fftfreq(np.size(x),d=dt)), np.fft.fft(x))
# ax[1].set_xlim(-1,1)
#
# # plt.axis('equal')
plt.show()
f.close()