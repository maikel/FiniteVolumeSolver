# /usr/bin/python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import os
import math

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

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
        data_array = np.array(file['data'][chunk[0]:chunk[1],:,:,:,:])
        time_array = np.array(file['times'][chunk[0]:chunk[1]])
    else:
        data_array = np.array(file['data'])
        time_array = np.array(file['times'])
    file.close()
    return (data_array, time_array)


def PlotFrame(fig, ax, cax, i0, i, data, func, times, dest):
    density = data[:,0,:,:,0]
    dims = [density.shape[2], density.shape[1]]
    mapped = np.reshape(func(data, i), dims)
    density_i = np.reshape(density[i,:,:], dims)
    mapped_i = np.where(density_i > 0, mapped, np.nan)

    mapped_f = np.nan_to_num(mapped_i, nan=np.finfo(float).max)
    vmin = 0.8e5

    mapped_f = np.nan_to_num(mapped_i, nan=np.finfo(float).min)
    vmax = 1.9e5

    xs = np.linspace(vmin, vmax, 4, endpoint=True)

    im = ax.imshow(mapped_i, vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, cax=cax, orientation='horizontal', ticks=xs)
    ax.set_title('Time = {:.6f}'.format(times[i]))
    os.makedirs(dest, exist_ok=True)
    filename = '{}/data_{:04d}.png'.format(dest, i0 + i)
    print('Save to {}.'.format(filename))
    fig.savefig(filename)

def ToPressure(data, i):
  return data[i, 16, :, :, 0]

def MakeImages(dest, source, steps=[], chunkSize=100):
    if steps:
        iStart = steps[0]
        iEnd = steps[1]
    else:
        iStart = 0
        iEnd = GetSize(source)
    for i in range(iStart, iEnd, chunkSize):
        thisChunkSize = min(chunkSize, iEnd - i)
        chunk = [i, i + thisChunkSize]
        data, times = Load(path, chunk)
        for j in range(0, thisChunkSize):
            fig, ax = plt.subplots()
            cax = fig.add_axes([0.27, 0.9, 0.5, 0.05])
            PlotFrame(fig, ax, cax, i, j, data, ToPressure, times, dest)
            cax.cla()
            ax.cla()
            fig.clf()

slurm_id = 'local'
path = '/srv/public/Maikel/FiniteVolumeSolver/build_3D-Release/MultiTube/Slices/Plenum.h5'
dest = '/srv/public/Maikel/PressureBlocked/{}'.format(slurm_id)
MakeImages(dest, path)