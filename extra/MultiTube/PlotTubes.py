import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import os
import math
from pathlib import Path


os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

def GetSize(path):
    file = h5py.File(path, mode='r')
    n = file["times"].shape[0]
    file.close()
    return n

def GetXLen(path):
    file = h5py.File(path, mode='r')
    n = file["data"].shape[2]
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

def Unzip(xs):
    return ([i for i, j in xs], [j for i, j in xs])

def LoadAll(paths, chunk=[]):
    return Unzip([Load(path, chunk) for path in paths])
    
def PlotTimeOverX(f, fig, axs, inputs, time):
    i = 0
    for (ax, input) in zip(axs, inputs):
        data = f(input)
        ax.imagesc(time, data, origin='lower')
        ax.set_title('Tube {}'.format(i))
        i = i + 1

def PressureField(input):
    return input[:, 14, :, 0, 0]

def H2Field(input):
    return input[:, 9, :, 0, 0]

def MakeImages(dest, sources, field_function = PressureField, steps=[], chunkSize=100, title=None, vmin=None, vmax=None):
    if steps:
        iStart = [steps[0] for source in sources]
        iEnd = [steps[1] for source in sources]
    else:
        iStart = [0 for source in sources]
        iEnd = [GetSize(source) for source in sources]
    n_cells = [GetXLen(source) for source in sources]
    print(n_cells)
    print(iEnd)
    all_fields = [np.zeros((i, n)) for i,n in zip(iEnd, n_cells)]
    all_times = [np.zeros(i) for i in iEnd]
    for i in range(iStart[0], iEnd[0], chunkSize):
        thisChunkSize = min(chunkSize, iEnd[0] - i)
        chunk = [i, i + thisChunkSize]
        datas, times = LoadAll(paths, chunk)
        fields = [field_function(data) for data in datas]
        for g, t, f, tt in zip(all_fields, all_times, fields, times):
            g[chunk[0]:chunk[1], :] = f
            t[chunk[0]:chunk[1]] = tt
    fig, axs = plt.subplots(2, 3, sharey='row', sharex='col', figsize=(10, 9))
    fig.tight_layout(pad=5.0)
    tube_id = 0
    for (ax, t, field) in zip(np.reshape(axs, 6), all_times, all_fields):
        xs = np.linspace(0.0, 1.5, field.shape[1])
        im = ax.imshow(field, origin='lower', interpolation='none', extent=(np.amin(xs), np.amax(xs), np.amin(t), np.amax(t)), aspect='auto', vmin=vmin, vmax=vmax)
        ax.set_title('Tube {}'.format(tube_id + 1))
        ax.set_xlabel('X-Position [m]')
        ax.set_ylabel('Time [s]')
        fig.colorbar(im, ax=ax)
        tube_id = tube_id + 1
    path = Path(dest)
    os.makedirs(path.parent, exist_ok=True)
    if title is not None:
        fig.suptitle(title, fontsize=16)
    fig.savefig(dest)
    fig.clf()

slurm_id = '5672798'
paths = ['/scratch/guttula/MultiTube_blocking/{}/MultiTube/Slices/Tube_{}.h5'.format(slurm_id, i) for i in range(0, 6)]
dest = '/home/guttula/FiniteVolumeSolver/extra/MultiTube/{}/H2.png'.format(slurm_id)
#dest = '/home/guttula/FiniteVolumeSolver/extra/MultiTube/{}/Pressure.png'.format(slurm_id)
#MakeImages(dest, paths, field_function=PressureField, chunkSize=1000, title='Pressure within the Tubes', vmin=1e4, vmax=2e5)
MakeImages(dest, paths, field_function=H2Field, chunkSize=1000, title='H2 within the Tubes')
