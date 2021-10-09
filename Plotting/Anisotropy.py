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
        energyArray.append(6.241509e21*1.0545718e-34*2*np.pi*np.abs(np.fft.fftfreq(np.size(x),d=dt)[np.argmax(np.abs(np.fft.fft(x)))]))
        i+=1
        path = basePath[:-4] + str(i) + basePath[-4:]
        print(f'file {i}')
    except IOError:
        if i==0:
           print('Error While Opening the file!')
        else:
           success = False





fig, ax = plt.subplots()
#


time = np.linspace(0,np.size(x)*dt,np.size(x))
#energy = data[3::N]

print(f'size: {np.size(x)}')
print(f'z difference: {np.max(z)-np.min(z)}')
print(np.size(np.fft.fft(x)))
print('plotting...')

# y = 9.274009994e-24*25/(1.38064852e-23*T)
# langevin = np.cosh(y)/np.sinh(y) - 1/y


print(f'precession energy: {4.135665538536e-12*np.abs(np.fft.fftfreq(np.size(x),d=dt)[np.argmax(np.abs(np.fft.fft(x)))])} MeV')

ax.plot(angleArray, 6.241509e21*2.002* 9.274009994e-24*10*np.cos(np.radians(angleArray)),'r')
ax.plot(angleArray, energyArray, 'bo')
ax.set_xlabel(r'$\theta$ [degrees]')
ax.set_ylabel(r'Energy [meV]')
ax.set_title("Anisotropy")
#
# # plt.axis('equal')
plt.show()
f.close()