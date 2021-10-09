import numpy as np
import cmath
import matplotlib.pyplot as plt

data = np.fromfile("C:/Users/timov/source/repos/CLaSSiC2.0/CLaSSiC2.0/data0.dat", dtype=np.double)

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

N = x.size
atoms = 12

x = x.reshape(atoms, N)
y = y.reshape(atoms, N)
z = z.reshape(atoms, N)



lattice_position = np.linspace(0,11,12)
print(lattice_position)
def scatter(q, lattice):
    t_0 = 1
    I_aa = np.zeros((3, int(N / 2)))
    ic_sum = np.zeros(3)
    q_dot_lattice = np.exp(1j * np.dot(q, lattice[0]))
    print(q_dot_lattice)
    for i in range(N):
        I_aa[0, :] += np.fft.fft(q_dot_lattice * x[i,:].T, n=N) \
                          .reshape((N,))[:int(N / 2)]
        I_aa[1, :] += np.fft.fft(q_dot_lattice * y[i,:].T, n=N) \
                          .reshape((N,))[:int(N / 2)]
        I_aa[2, :] += np.fft.fft(q_dot_lattice * z[i,:].T, n=N) \
                          .reshape((N,))[:int(N / 2)]

        ic_sum[0] += cmath.exp(-1j * np.dot(q, lattice[i])) * x[i, 0]
        ic_sum[1] += cmath.exp(-1j * np.dot(q, lattice[i])) * y[i, 0]
        ic_sum[2] += cmath.exp(-1j * np.dot(q, lattice[i])) * z[i, 0]

    for i in range(3):
        I_aa[i, :] *= ic_sum[i]
        I_aa[i, :] = np.square(np.abs(I_aa[i, :]), 2)

    frequencies = 2*np.pi*np.fft.fftfreq(N, dt)[:int(N / 2)]
    energies = 4.1357e-15 * frequencies



scatter([0.2,0,0],[[1,0,1],[-1,0,1]])