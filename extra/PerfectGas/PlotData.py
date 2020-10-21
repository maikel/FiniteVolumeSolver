import h5py
import numpy
import matplotlib.pyplot as plt

def LoadDataAndTimes(path):
  file = h5py.File(path, 'r')
  data = numpy.array(file['data'])
  times = numpy.array(file['times'])
  file.close()
  return data, times

def ImageShow(ax, field, extent):
  im=ax.imshow(field, origin='lower', aspect='auto', interpolation='none', extent=extent)
  # plt.colorbar(im,ax=ax)

data1,times1 = LoadDataAndTimes('WithDiffusion.h5')
data2,times2 = LoadDataAndTimes('NoDiffusion.h5')
fig, axs = plt.subplots(1,2)

rho1 = data1[:, 0, :, 0]
rhou1 = data1[:, 1, :, 0]
rhoE1 = data1[:, 2, :, 0]
p1 = data1[:, 3, :, 0]

rho2 = data2[:, 0, :, 0]
rhou2 = data2[:, 1, :, 0]
rhoE2 = data2[:, 2, :, 0]
p2 = data2[:, 3, :, 0]

ImageShow(axs[0], rho1, (0.0, 2.0, times1[0], times1[-1]))
axs[0].set_title('Density With Diffusion')

ImageShow(axs[1], rho2, (0.0, 2.0, times2[0], times2[-1]))
axs[1].set_title('Density Without Diffusion')

fig.savefig('Density.png')