import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
# print(FVS_path) 
sys.path.append(FVS_path+'/extra/')
# import amrex.plotfiles as da
from amrex.plotfiles import h5_load_timeseries

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


sys.path.append(inputFilePath)
from SEC_Plenum_Arrhenius import t_ref, ControlOptions, T_ref

controlState = "{}/ControlState.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization'.format(outPath)



os.makedirs(output_path, exist_ok=True)


datas, times, datas_dict = h5_load_timeseries(controlState)

def getSubKeyList(substring):
   return [ key for key in datas_dict.keys() if substring in key ]
   
f, axs = plt.subplots(nrows=3, ncols=2, figsize=(17/2,20/2) )
f.suptitle('Control Station')

plotKeyList = ['mass_flow', 'power', 'pressure', 'temperature', 'rpm', 'fuel_consumption_rate']

datas[:, datas_dict['current_rpm']] = datas[:, datas_dict['current_rpm']] / ControlOptions['rpmmax']


data_units = ['', '', ' [bar]', ' [K]', ' [-]', '']

for i, ax, subKey in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList):
   subKeyList = getSubKeyList(subKey)
   for key in subKeyList:
      lab = key.replace(subKey, '')
      lab = lab.replace('_', ' ')
      if (i==0 and 'comp' in lab):
         lab = lab.replace('compressor', 'comp. plenum')

      if 'temperature' in subKey:
         datas[:, datas_dict[key]] *=T_ref
      
      ax.plot( times, datas[:, datas_dict[key]], label=lab )
      
      if i>3:
         pass
      else:
         ax.legend()
      ylab = subKey.replace('_', ' ')
      if 'rpm' in ylab:
         ylab='rpm / rpmmax'
      ylab += data_units[i]
      ax.set(xlabel='time', ylabel=ylab, xlim=(0,None))

f.subplots_adjust(hspace=0.3, wspace=0.3)
f.savefig('{}/control_state.png'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)
