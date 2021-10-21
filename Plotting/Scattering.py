import helper
import numpy as np
import cmath
import matplotlib.pyplot as plt

param, x, y, z = helper.getData()
fNum = 0

frequencies = np.fft.fftfreq(param[fNum]["steps"], param[fNum]["dt"])
energies = 4.1357e-12 * frequencies
# subN = int(N/2/100)
freqN = 12
maxEnergyIndex = np.argmax(energies> 6)
print(maxEnergyIndex)
lattice_position = np.linspace(0,11e-10,12)
lattice_position = []
fourierLength = int(param[fNum]["steps"]/2)
# fourierLength = maxEnergyIndex

for i in range(12):
    lattice_position.append([i,0,0])
lattice_position = np.array(lattice_position).reshape(12,3)


def scatter(q, position):
    I_aa = np.zeros((3, fourierLength), dtype=complex)
    ic_sum = np.zeros(3, dtype=complex)
    for i in range(param[fNum]["atoms"]):
        q_dot_lattice = np.exp(1j * np.dot(q, position[i]))
    
        I_aa[0, :] += np.fft.fft(q_dot_lattice * x[fNum, i, :].T)[:fourierLength]
        I_aa[1, :] += np.fft.fft(q_dot_lattice * y[fNum, i, :].T)[:fourierLength]
        I_aa[2, :] += np.fft.fft(q_dot_lattice * z[fNum, i, :].T)[:fourierLength]

        ic_sum[0] += cmath.exp(-1j * np.dot(q, position[i])) * x[fNum, i, 0]
        ic_sum[1] += cmath.exp(-1j * np.dot(q, position[i])) * y[fNum, i, 0]
        ic_sum[2] += cmath.exp(-1j * np.dot(q, position[i])) * z[fNum, i, 0]

    for i in range(3):
        I_aa[i, :] *= ic_sum[i]
        I_aa[i, :] = np.power(np.abs(I_aa[i, :]), 2)



    return I_aa.real


I_total = np.zeros((freqN, maxEnergyIndex))
# subEnergies = np.zeros(100)
# for i in range(100-1):
#     subEnergies[i] = np.sum((energies[:maxEnergyIndex])[i*subN:(i+1)*subN])/subN
# print(subEnergies)

q = np.array([0,0,0])
for i in range(freqN):
    print(f'q {i}')
    q = np.vstack((q, np.array([np.pi*(2*i/freqN-1),0,0])))
    I_aa = scatter(q[-1,:],lattice_position)

    print(energies[np.argmax(I_aa[0,:])])
    print(np.max(I_aa[0, :]))
    # plt.plot(energies[:fourierLength], I_aa[1,:])
    
    # for j in range(100-1):
        # I_total[i, j] = np.sum(np.flip(I_aa[1,:maxEnergyIndex])[j*subN:(j+1)*subN].real)/subN
    I_total[i, :] = I_aa[0,:maxEnergyIndex]/np.max(I_aa[0,:maxEnergyIndex])

    # plt.plot(np.flip(energies[:maxEnergyIndex]), I_total[i, :], '.')


# I_total = I_total**10
# I_total /= np.max(I_total)
print(f'max: {np.max(I_total)}')
print(f'shape: {I_total.shape}')



# plt.hist2d(I_total[:,0],I_total[], bins=[20, 100], range=[[0,20*0.25],[np.min(energies),np.max(energies)]])
# Construct the plot
# plt.figure(figsize=(10, 9), tight_layout=True)

plt.plot(q[1:,0], helper.constants["J_to_meV"] * 4 * param[fNum]["J"]*3.5*(1-np.cos(q[1:,0])), 'r')
for i in range(I_total.shape[0]):
    if np.max(I_total[i, :]) > 0.2:
        plt.plot(q[1+i,0], energies[np.argmax(I_total[i,:])], 'b.')

# extent = [np.min(q)-0.5*1/freqN, np.max(q), 0, 6]
# print(extent)
# # print(I_total[:,:maxEnergyIndex].shape)
# im = plt.imshow(I_total.T, extent=extent, cmap='Blues')
# ax = plt.gca()
# ax.set_aspect(abs(extent[1] - extent[0]) / abs(extent[3] - extent[2]))
# cbar = ax.figure.colorbar(im, ax=ax)

# im.set_clim(vmin=0.0, vmax=1.0)

# print(	1.602e22 * 4 * J * 3.5 * (1-np.cos(q[1:,0])))

# Prettify
# Start with the labels
# fontsize = 18
# q_labels = ['0', '0', '0']
# q_labels[direction_number] = 'X'
# plt.xlabel(f'$q = ({q_labels[0]}, {q_labels[1]}, {q_labels[2]}) [Å^{"{-1}"}]$', fontsize=fontsize)
# plt.ylabel('Energy [meV]', fontsize=fontsize)
# cbar.set_label('Normalized intensity [A.U.]', fontsize=fontsize)

# # Now the ticks
# plt.yticks(fontsize=16)
# plt.xticks(fontsize=16)
# cbar.set_ticks([])
plt.show()