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


# import importlib
# inputFile = importlib.import_module('SEC_Plenum')
# print(inputFile.y0s)
# print(inputFile.Area)
sys.path.append(inputFilePath)
from SEC_Plenum_Arrhenius import y0s, Area, tube_n_cells, p_ref, T_ref, rho_ref, Output
# print(y0s)
# print(Area)


controlState = "{}/ControlState.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization'.format(outPath)

def h5_load(path, variables):
   file = h5py.File(path, mode='r')
   strings = list(file['fields'].asstr())
   indices = [strings.index(var) for var in variables]
   dictionary = dict( zip(strings, indices) )
   shape = file['data'].shape
   datas = np.zeros((shape))
   for i ,var in enumerate(indices):
      datas[:,i] = file['data'][:, var]
   time = np.array(file['times'])
   file.close()
   return datas, time, dictionary


# Get nsteps from controlState
file = h5py.File(controlState, mode='r')
nsteps = file['data'].shape[0]
file.close()


def PrintProgress(i):
  progress = int(100.0 * float(i) / (nsteps - 1))
  print('[{:3d}%] Reading slice [{}/{}]'.format(progress, i + 1, nsteps))

os.makedirs(output_path, exist_ok=True)


variables = ["compressor_mass_flow", "compressor_power", "compressor_pressure",
   "compressor_temperature", "current_rpm", "fuel_consumption",
   "power_out", "turbine_mass_flow", "turbine_power",
   "turbine_pressure", "turbine_temperature"]

datas, times, datas_dict = h5_load(controlState, variables)
# print(datas.shape)
# print(datas[:, datas_dict['current_rpm']])

   
f, axs = plt.subplots(nrows=3, ncols=2, figsize=(17/2,20/2) )
f.suptitle('control_state')

tmpList = ['compressor_mass_flow', 'turbine_mass_flow']
for el in tmpList:
   axs[0,0].plot( times, datas[:, datas_dict[el]], label=el )
   axs[0,0].legend()

tmpList = ['compressor_power', 'turbine_power']
for el in tmpList:
   axs[0,1].plot( times, datas[:, datas_dict[el]], label=el )
   axs[0,1].legend()

tmpList = ['compressor_pressure', 'turbine_pressure']
for el in tmpList:
   axs[1,0].plot( times, datas[:, datas_dict[el]], label=el )
   axs[1,0].legend()

tmpList = ['compressor_temperature', 'turbine_temperature']
for el in tmpList:
   axs[1,1].plot( times, datas[:, datas_dict[el]], label=el )
   axs[1,1].legend()

tmpList = ['current_rpm']
for el in tmpList:
   axs[2,0].plot( times, datas[:, datas_dict[el]], label=el )
   axs[2,0].legend()

tmpList = ['fuel_consumption']
for el in tmpList:
   axs[2,1].plot( times, datas[:, datas_dict[el]], label=el )
   axs[2,1].legend()

for ax in axs.flatten():
   ax.set(xlabel='time')

f.savefig('{}/control_state.png'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)
