import matplotlib.pyplot as plt
import os
import h5py
import numpy as np

os.environ['HDF5_USE_FILE_LOCKING'] = 'TRUE'

def GetSize(path):
    file = h5py.File(path, mode='r')
    n = file["times"].shape[0]
    file.close()
    return n

def Load(path, chunk=[]):
    # @brief Reads grid data and time data from a specified hdf5 file.
    # @param path The path of the HDF5 file
    # @param dataset A string name of the main dataset containing the states data
    # @param times A string name of the time dataset containing all time points
    # @param chunk An optional [start, end] list that specifies how many time steps to load
    file = h5py.File(path, "r")
    if chunk:
        data_array = np.array(file['data'][chunk[0]:chunk[1],0,:,:])
        time_array = np.array(file['times'][chunk[0]:chunk[1]])
    else:
        data_array = np.array(file['data'][:, 0, :, :])
        time_array = np.array(file['times'])
    file.close()
    return data_array #, time_array)

path = '/home/amrex/Divider_c24_800.h5'
output_path = '/home/amrex/Figures/Hdf5'
os.makedirs(output_path, exist_ok=True)

file = h5py.File(path, mode='r')
rho_data = np.array(file['data'][:, 0, :, :])
rhou_data = np.array(file['data'][:, 1, :, :])
rhov_data = np.array(file['data'][:, 2, :, :])
rhoE_data = np.array(file['data'][:, 3, :, :])
p_data = np.array(file['data'][:, 4, :, :])
c_data = np.array(file['data'][:, 5, :, :])
alpha = np.array(file['data'][:, -1, :, :])
times = np.array(file['times'])
dx = file['data'].attrs['cell_size']
xlower = file['data'].attrs['xlower']
file.close()

nx = [rho_data.shape[1], rho_data.shape[2]]
xupper = xlower + dx * nx

vmin = None #0.4
vmax = None #2.0

#ns = range(times.shape[0]) # [i for i in range(data.shape[0]) if i % 5 == 0]
ns = [i for i in range(times.shape[0]) if i % 5 == 0]
for i, n in enumerate(ns):
    p = p_data[n,:,:]
    rho = rho_data[n,:,:]
    extent = np.array([xlower[0], xupper[0], xlower[1], xupper[1]])
    volfrac = alpha[n, :, :]
    p = np.reshape(p, [nx[1], nx[0]])
    rho = np.reshape(rho, [nx[1], nx[0]])
    volfrac = np.reshape(volfrac, [nx[1], nx[0]])
    p = np.where(volfrac == 0.0, np.nan, p)
    rho = np.where(volfrac == 0.0, np.nan, rho)
    f, axs = plt.subplots(nrows=2, ncols=1, figsize=(25, 20))
    f.suptitle('Time = {:.2e}'.format(times[n]))
    im_p = axs[0].imshow(p, origin='lower', vmin=0.9e5, vmax=1.4e5, interpolation='none', extent=extent, cmap='jet')
    axs[0].set_title('Pressure')
    plt.colorbar(im_p, ax=axs[0])
    # Temperature image
    Rspec = 287.058
    T = p / rho / Rspec
    im_T = axs[1].imshow(T, origin='lower', vmin=290.0, vmax=315.0, interpolation='none', extent=extent, cmap='jet')
    axs[1].set_title('Temperature')
    plt.colorbar(im_T, ax=axs[1])
    # save figure
    plt.tight_layout()
    figure_path = '{}/Figure{:04d}.png'.format(output_path, i)
    print('Save figure to {}.'.format(figure_path))
    f.savefig(figure_path)
    f.clear()
    plt.close(f)
