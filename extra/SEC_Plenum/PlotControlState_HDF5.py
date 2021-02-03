import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
# print(FVS_path) 
sys.path.append(FVS_path+'/extra/')

import yt
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import h5py
import itertools

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

dataPath = FVS_path+"/build_2D-Release/average_massflow"
inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"


# sys.path.append(inputFilePath)
# from SEC_Plenum_Arrhenius import y0s, Area, tube_n_cells, p_ref, T_ref, rho_ref, Output



controlState = "{}/ControlState.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization'.format(outPath)

def h5_load_timeseries(path):
   file = h5py.File(path, mode='r')
   strings = list(file['fields'].asstr())
   indices = [strings.index(var) for var in strings]
   dictionary = dict( zip(strings, indices) )
   shape = file['data'].shape
   datas = np.zeros((shape))
   for i ,var in enumerate(indices):
      datas[:,i] = file['data'][:, var]
   time = np.array(file['times'])
   file.close()
   return datas, time, dictionary


os.makedirs(output_path, exist_ok=True)


datas, times, datas_dict = h5_load_timeseries(controlState)

def getSubKeyList(substring):
   return [ key for key in datas_dict.keys() if substring in key ]
   
f, axs = plt.subplots(nrows=3, ncols=2, figsize=(17/2,20/2) )
f.suptitle('control_state')

plotKeyList = ['mass_flow', 'power', 'pressure', 'temperature', 'rpm', 'fuel']

for ax, subKey in zip(axs.flatten(), plotKeyList):
   subKeyList = getSubKeyList(subKey)
   for key in subKeyList:
      ax.plot( times, datas[:, datas_dict[key]], label=key )
      ax.legend()
      ax.set(xlabel='time')

f.subplots_adjust(hspace=0.3)
f.savefig('{}/control_state.png'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)
