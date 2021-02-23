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

# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1])
   inputFilePath = dataPath
else:
   dataPath = FVS_path+"/build_2D-Release/average_massflow"
   inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"

# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = True

sys.path.append(inputFilePath)
from SEC_Plenum_Arrhenius import t_ref, ControlOptions, T_ref

outPath = dataPath
output_path = '{}/Visualization'.format(outPath)

plt.style.use('seaborn')
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "font.size": 9,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7
}
plt.rcParams.update(tex_fonts)

os.makedirs(output_path, exist_ok=True)

filename_basic = "{}/ControlState.h5".format(dataPath)

if RESTARTEDSIMULATION:
   import glob
   fileNameList = [filename_basic]

   ###### collect data begin
   # check if other h5.* files exist and append them to the list
   otherFiles = glob.glob("{}.*".format(filename_basic))
   if otherFiles:
      fileNameList.append( *otherFiles )

   print(fileNameList)

   # Read in data
   # Attention the last file is the latest!!
   # for example we have Filename.h5 and Filename.h5.1 the last one contains the first data!
   print("Read in data from {}".format(fileNameList[-1]))
   datas, times, datas_dict = h5_load_timeseries(fileNameList[-1])

   for filename in reversed(fileNameList[:-1]):
      print("Read in data from {}".format(filename))
      data, time, _ = h5_load_timeseries(fileNameList[0])
      datas = np.concatenate((datas, data))
      times = np.concatenate((times, time))
else:
   print("Read in data from {}".format(filename_basic))
   datas, times, datas_dict = h5_load_timeseries(filename_basic)


def getSubKeyList(substring):
   return [ key for key in datas_dict.keys() if substring in key ]
   
f, axs = plt.subplots(nrows=4, ncols=2, figsize=(23/2,20/2) )
f.suptitle('Control Station')


plotKeyList = ['mass_flow', 'power', 'pressure', 'temperature', 'rpm', 'fuel_consumption_rate', 'efficiency', 'SEC_Mode']
# print(datas_dict)
datas[:, datas_dict['current_rpm']] = datas[:, datas_dict['current_rpm']] / ControlOptions['rpmmax']
print( "maximum compressor_pressure = {}".format(np.max(datas[:, datas_dict['compressor_pressure']])) )

# calculate the carnot efficiency following Kleins SEC_Code
eff_carnot = 1.0 - 1.0/datas[:, datas_dict['compressor_temperature']]

data_units = ['', '', ' [bar]', ' [K]', ' [-]', '', '', '']

for i, ax, subKey in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList):
   subKeyList = getSubKeyList(subKey)
   for key in subKeyList:
      if i>5: # don't destroy the labels from the last two plots
         lab = key
      else:
         lab = key.replace(subKey, '')
         lab = lab.replace('_', ' ')
      
      if (i==0):
         if ('comp' in lab):
            lab = lab.replace('compressor', 'comp. plenum')
         elif ('turb' in lab):
            lab = lab.replace('turbine', 'turb. plenum')
      
      if 'temperature' in subKey:
         datas[:, datas_dict[key]] *=T_ref
      
      ax.plot( times, datas[:, datas_dict[key]], label=lab )
      
      if i>3:
         pass
      else:
         ax.legend()
      
      if 'efficiency' in subKey:
         ax.plot(times, eff_carnot, label='Carnot efficiency')
         ax.legend()
         ax.set(ylim=(0,0.5))
      
      ylab = subKey.replace('_', ' ')
      if 'rpm' in ylab:
         ylab='rpm / rpmmax'
      ylab += data_units[i]
         
      ax.set(xlabel='time', ylabel=ylab, xlim=(times[0], None))
      ax.grid(True)

f.subplots_adjust(hspace=0.3, wspace=0.3)
f.savefig('{}/control_state.png'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)