import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.plotfiles as da

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import h5py

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1])
   inputFilePath = dataPath
else:
   dataPath = FVS_path+"/build_2D-Release/average_massflow"
   inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"

try:
  inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
  inputfileName = 'SEC_Plenum_Arrhenius.py'

da.import_file_as_module(inputFilePath+inputfileName, 'inputfile')
from inputfile import Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref
from inputfile import D as diameter_tube

try:
  from inputfile import T_ref, n_tubes  
except:
  from inputfile import T_ref
  n_tubes=1

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
plt.rcParams.update({'axes.grid' : False})

#--------------------------------------------------------------------------

plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization/PVDiagramm/'.format(outPath)

# Get nsteps from plenum
file = h5py.File(plenum, mode='r')
nsteps = file['data'].shape[0]
file.close()



def PrintProgress(i):
  progress = int(100.0 * float(i) / (nsteps - 1))
  print('[{:3d}%] Reading slice [{}/{}]'.format(progress, i + 1, nsteps))

os.makedirs(output_path, exist_ok=True)

tube_paths = ['{}/Tube{}.h5'.format(dataPath, tube) for tube in range(n_tubes)]

def stackTubeDataTo2D(Tube_datalist):
   # all Tubedata is 1D but for contourf we need at least 2D data. so simply stack twice the 1d array
   for i, el in enumerate(Tube_datalist):
      el = np.squeeze(el)
      Tube_datalist[i] = np.stack((el,el))
   return Tube_datalist


def find_nearest(array, value, oldID):
  array = np.asarray(array)
  return (np.abs(array - value)).argmin()
  # while True:
  #   idx = (np.abs(array - value)).argmin()
  #   if (idx==(array.shape[-1]-1)) or (idx>=oldID):
  #     break
  #   else:
  #     print(idx, array[idx], value, oldID)
  #     array = np.delete(array, idx)
  #     oldID-=1

  # return idx

#---------------------------
test_skalar = 610.0 # initial passive scalar value 
tplotmin = 200.0
tplotmax = 400.0
# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False
#-------------------------------------------

## list to collect all the data
## [0] --> pressure
## [1] --> density (later specific volume)
passive_skalar = [[],[]]

for tube_id in range(n_tubes):
  print("[Tube{}] Plotting Tube data".format(tube_id))
  filename_basic = dataPath+'/Tube{}.h5'.format(tube_id)
  # datas, times, datas_dict = da.h5_load_timeseries(filename_basic)
  extent_1d = da.h5_load_get_extent_1D(filename_basic)
  
  if RESTARTEDSIMULATION:
   import glob
   fileNameList = [filename_basic]

   ###### collect data begin
   # check if other h5.* files exist and append them to the list
   otherFiles = glob.glob("{}.*".format(filename_basic))
   if otherFiles:
      fileNameList.append( *otherFiles )

   print("[Tube{}] fileNameList = {}".format(tube_id, fileNameList))

   # Read in data
   # Attention the last file is the latest!!
   # for example we have Filename.h5 and Filename.h5.1 the last one contains the first data!
   print("[Tube{}] Read in data from {}".format(tube_id, fileNameList[-1]))
   datas, times, datas_dict = da.h5_load_timeseries(fileNameList[-1])

   for filename in reversed(fileNameList[:-1]):
      print("[Tube{}] Read in data from {}".format(tube_id, filename))
      data, time, _ = da.h5_load_timeseries(fileNameList[0])
      datas = np.concatenate((datas, data))
      times = np.concatenate((times, time))
  else:
    print("[Tube{}] Read in data from {}".format(tube_id, filename_basic))
    datas, times, datas_dict = da.h5_load_timeseries(filename_basic)


  datas = np.squeeze(datas) # remove last axis
  print("[Tube{}] data shape from tube is {} = (NTimePoints, NVariables, NCells)".format(tube_id, datas.shape))
  # print(datas_dict)

  # optional slicing in time-dimension
  t_index_array = (times>=tplotmin) & (times<=tplotmax)
 
  # check if index array is empty
  if not np.any(t_index_array):
    raise IndexError('time index array is empty!')

  rho_data = datas[t_index_array, datas_dict['Density'], :]
#   rhou_data = datas[t_index_array, datas_dict['Momentum'], :]
#   rhoE_data = datas[t_index_array, datas_dict['Energy'], :]
#   rhoF_data = datas[t_index_array, datas_dict['Species'], :]
  p_data = datas[t_index_array, datas_dict['Pressure'], :]
#   c_data = datas[t_index_array, datas_dict['SpeedOfSound'], :]
#   T_data = p_data / rho_data
#   F_data = rhoF_data / rho_data
#   Ma = rhou_data / rho_data / c_data
  
  if not 'PassiveScalars' in datas_dict:
   raise KeyError('No PassiveScalars data could be found!')
   
  rhoX_data = datas[t_index_array, datas_dict['PassiveScalars'], :]
  skalarX_data = rhoX_data / rho_data
  indexed_time = times[t_index_array]
  
  def valueCheck(value, array):
    if (value < np.min(array)):
      print('scalar value: {} is lower than min scalar value: {}'.format(value, np.min(array)))
      return False
      # raise ValueError('scalar value: {} is lower than min scalar value: {}'.format(value, np.min(array)))
    elif (value > np.max(array)):
      print('scalar value: {} is greater than max scalar value: {}'.format(value, np.max(array)))
      return False
      # raise ValueError('scalar value: {} is greater than max scalar value: {}'.format(value, np.max(array)))
    else:
      return True

  if not valueCheck(test_skalar, skalarX_data):
    raise ValueError("value Check failed!")
  
  # if (test_skalar > np.max(skalarX_data)):
  #   raise ValueError('scalar value: {} is greater than max scalar value: {}'.format(test_skalar, np.max(skalarX_data)))
  # elif (test_skalar < np.min(skalarX_data)):
  #   raise ValueError('scalar value: {} is lower than min scalar value: {}'.format(test_skalar, np.min(skalarX_data)))

  print(skalarX_data.shape)
  print(np.min(skalarX_data), np.max(skalarX_data))

  test = [[],[], []]

  for step in range(skalarX_data.shape[0]):
    print("Search for {}".format(test_skalar))
    current_time = indexed_time[step]
    if not valueCheck(test_skalar, skalarX_data[step,:]):
      continue
    
    if not test[2]:
      new_idx = find_nearest(skalarX_data[step,:], test_skalar, 0)
    else:
      new_idx = find_nearest(skalarX_data[step,:], test_skalar, test[2][-1])
    print("found index {}".format(new_idx))
    print("current time is: {}".format(current_time))
   
    test_skalar = skalarX_data[step,new_idx]
    print("value from index is {}".format(test_skalar))

    passive_skalar[0].append(p_data[step,new_idx])
    passive_skalar[1].append(rho_data[step,new_idx])
    test[0].append(test_skalar)
    test[1].append(current_time)
    test[2].append(new_idx)

    if new_idx==skalarX_data.shape[1]-1:
      print("scalar has left the tube!")
      print("current time = {}".format(current_time))
      break
  passive_skalar = np.asarray_chkfinite(passive_skalar)
  # print(passive_skalar)

  f, ax = plt.subplots(nrows=1, ncols=1)
  ax.plot(passive_skalar[1]**(-1), passive_skalar[0], '-x')
  ax.set(xlabel="specific volume", ylabel="pressure")
  ax.set(title="scalar value is approx: {}".format(np.mean(test[0])))
  f.savefig('{}/pvDiagramm_tubeID{}.png'.format(output_path, tube_id), bbox_inches='tight')

  f, ax = plt.subplots(nrows=1, ncols=1)
  ax.plot(test[1], test[0], 'x')
  ax.set(ylabel="scalar value", xlabel="time")
  f.savefig('{}/pvScalarVsTime_tubeID{}.png'.format(output_path, tube_id), bbox_inches='tight')
  #
  f, ax = plt.subplots(nrows=1, ncols=1)
  ax.plot(test[1], test[2], 'x')
  ax.set(ylabel="scalar index", xlabel="time")
  f.savefig('{}/pvIDVsTime_tubeID{}.png'.format(output_path, tube_id), bbox_inches='tight')


