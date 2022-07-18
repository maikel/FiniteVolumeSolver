# /usr/bin/python

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
    file = h5py.File(path, "r", swmr=True)
    if chunk:
        data_array = np.array(file['data'][:,:,chunk[0]:chunk[1]])
        time_array = np.array(file['times'][chunk[0]:chunk[1]])
    else:
        data_array = np.array(file['data'])
        time_array = np.array(file['times'])
    file.close()
    return (data_array, time_array)

def FindLargestGradient(x):
  dx = x[1:-1] - x[0:-2]
  return np.argmax(dx > 1000)

def MakePlot(dest, path, tube_id=0, variable=0, steps=[], time_interval=[]):
    data, times = Load(path, steps)
    if time_interval:
      iStart = np.argmax(times > time_interval[0])
      rtimes = times[::-1]
      iEnd = len(times) - np.argmax(rtimes < time_interval[1]) - 1
      data = data[:,:,iStart:iEnd]
      times = times[iStart:iEnd]
    fig, ax = plt.subplots()
    ax.plot(times, data[variable, 0 + tube_id*5, :])
    ax.plot(times, data[variable, 1 + tube_id*5, :])
    ax.plot(times, data[variable, 2 + tube_id*5, :])
    ax.plot(times, data[variable, 3 + tube_id*5, :])
    ax.plot(times, data[variable, 4 + tube_id*5, :])
    ax.set(xlim=[34.8e-3, 40e-3])
    fig.savefig(dest)
    fig.clf()
    return data, times

slurm_id = '5343973'
# source = '/scratch/guttula/MultiTube_blocking/{}/MultiTube/Probes/Plenum.h5'.format(slurm_id)
source = '/srv/public/Maikel/FiniteVolumeSolver/build_3D-Release/MultiTube2/Probes/Plenum.h5'
#source = "/scratch/guttula/MultiTube/5097856/MultiTube/Probes/Plenum.h5"
dest = "/srv/public/Maikel/FiniteVolumeSolver/extra/MultiTube/Pressure.png"
variable = 16
tube_id = 0

data, times = MakePlot(dest, source, tube_id=tube_id, variable=variable)

# imax = [FindLargestGradient(data[16,i,:]) for i in range(0, 4)]
# print('Largest gradient at i={} and t={}s'.format(imax[0], times[imax[0]]))
# times_normalized = times - times[imax[0]]

# fig, ax = plt.subplots()
# ax.plot(times_normalized, data[variable, 0 + tube_id*5, :])
# ax.plot(times_normalized, data[variable, 1 + tube_id*5, :])
# ax.plot(times_normalized, data[variable, 2 + tube_id*5, :])
# ax.plot(times_normalized, data[variable, 3 + tube_id*5, :])
# ax.plot(times_normalized, data[variable, 4 + tube_id*5, :])
# ticks = [times_normalized[i] for i in imax]
# ax.set_xticks(ticks[1:])
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Static Pressure [Pa]')
# ax.set_title('Single Tube Case: the tube is filled with fuel up to 90%')
# ax.grid()

# dest = "/home/guttula/FiniteVolumeSolver/extra/MultiTube/FillRatio_90/pressure_90.png"
# fig.savefig(dest)
# fig.clf()