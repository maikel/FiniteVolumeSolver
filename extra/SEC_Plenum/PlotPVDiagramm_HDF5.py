import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
from amrex.other import import_file_as_module

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import h5py

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# check cli
if len(sys.argv)<2:
   errMsg = ('Not enough input arguments!\n'
               +'\tfirst argument must be dataPath!')
   raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

try:
  inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
  inputfileName = 'SEC_Plenum_Arrhenius.py'

import_file_as_module(os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref, L_ref, R_ref, R, gamma
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
plt.rcParams.update({'axes.grid' : True})

#--------------------------------------------------------------------------
# dicts and function
PerfectGasConstants = {
  'Rspec': R, # [-]
  'gamma': gamma, # [-]
}

ParameterNonDim = {
  'time' : t_ref, # [s]
  'length': L_ref, # [m]
  'pressure' : p_ref, # [Pa]
  'density' : rho_ref, # [kg/m^3]
  'temperature': T_ref, # [K]
  'velocity': u_ref, # [m/s]
  'SpecGasConstant': R_ref, # [J/(kg*K)]
}

def getSpecGasConstant(Dimless=True):
  Rspec = PerfectGasConstants['Rspec']
  if not Dimless:
    Rspec *= ParameterNonDim['SpecGasConstant']
  return Rspec

def getHeatCapacityConstantVolume(Dimless=True):
  gammaMinusOne = PerfectGasConstants['gamma']-1.
  Rspec = getSpecGasConstant(Dimless)
  return Rspec/gammaMinusOne

def getHeatCapacityConstantPressure(Dimless=True):
  gamma = PerfectGasConstants['gamma']
  gammaMinusOne = gamma-1.
  Rspec = getSpecGasConstant(Dimless)
  return Rspec*gamma/gammaMinusOne

def findNearest1D(array, value):
  array = np.asarray(array)
  sortedIndices = np.argsort(np.abs(array - value))
  
  idPos = []
  idNeg = []
  # greater than value
  for id in sortedIndices:
    if array[id]>=value:
      idPos.append(id)
      break
  # less than value
  for id in sortedIndices:
    if array[id]<value:
      idNeg.append(id)
      break
  if (not idPos):
    raise AssertionError(f"no value is greater than {value}!")
  elif (not idNeg):
    raise AssertionError(f"no value is less than {value}!")
  
  # get rid of list type
  idPos, = idPos
  idNeg, = idNeg

  return idNeg, idPos
  # return (np.abs(array - value)).argmin()

def computeEntropie(T, rho, s0=0., Dimless=True):
  cv = getHeatCapacityConstantVolume(Dimless)
  Rspec = getSpecGasConstant(Dimless)
  print(Dimless, Rspec, cv)
  s = np.zeros_like(rho)
  s[0] = s0 # set first element
  sdiff = (cv * np.log(T[1:]/T[:-1]) 
            + Rspec * np.log(rho[:-1]/rho[1:])
          )
  for i in range(s.shape[0]-1):
    s[i+1] = s[i]+sdiff[i]
  return s

def computeTemperature(p, rho, Dimless=True):
  Rspec = getSpecGasConstant(Dimless)
  return p/(rho*Rspec)

def valueCheck(value, array):
    if (value < np.min(array)):
      print('scalar value: {} is lower than min scalar value: {}'.format(value, np.min(array)))
      return False
    elif (value > np.max(array)):
      print('scalar value: {} is greater than max scalar value: {}'.format(value, np.max(array)))
      return False
    else:
      return True

def interpolate1D(value, x, y, ids):
      data = np.interp( value, 
                [x[ids[0]], x[ids[1]]], 
                [y[ids[0]], y[ids[1]]] 
               )
      return data

def simplePlot(strX, strY, path, markerstyle='x', tube_id=0):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(valueDict[strX]['data'], valueDict[strY]['data'], 
              '-', marker=markerstyle)
    ax.set(
      xlabel=f"{valueDict[strX]['label']} {valueDict[strX]['dim']}", 
      ylabel=f"{valueDict[strY]['label']} {valueDict[strY]['dim']}"
      )
    fig.savefig('{}/{}{}Diagramm_tubeID{}.png'.format(path, valueDict[strY]['symbol'], valueDict[strX]['symbol'], tube_id), bbox_inches='tight')
#-------------------------------------------------

plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization/PVDiagramm/'.format(outPath)

os.makedirs(output_path, exist_ok=True)

tube_paths = ['{}/Tube{}.h5'.format(dataPath, tube) for tube in range(n_tubes)]

#-------------------------------------------

#---------------------------
test_scalar = 820.0 # initial passive scalar value 610, 805
# best: 815, 820, 
# super: 840 
tplotmin = 350. #200.0
tplotmax = 400.0
Dimless = True
# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False
#-------------------------------------------

## list to collect all the data
## [0] --> pressure
## [1] --> density (later specific volume)
passive_scalar = [[],[]]

for tube_id in range(n_tubes):
  print("[Tube{}] Plotting Tube data".format(tube_id))
  filename_basic = dataPath+'/Tube{}.h5'.format(tube_id)
  # datas, times, datas_dict = io.h5_load_timeseries(filename_basic)
  extent_1d = io.h5_load_get_extent_1D(filename_basic)
  
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
   datas, times, datas_dict = io.h5_load_timeseries(fileNameList[-1])

   for filename in reversed(fileNameList[:-1]):
      print("[Tube{}] Read in data from {}".format(tube_id, filename))
      data, time, _ = io.h5_load_timeseries(fileNameList[0])
      datas = np.concatenate((datas, data))
      times = np.concatenate((times, time))
  else:
    print("[Tube{}] Read in data from {}".format(tube_id, filename_basic))
    datas, times, datas_dict = io.h5_load_timeseries(filename_basic)


  datas = np.squeeze(datas) # remove last axis
  print("[Tube{}] data shape from tube is {} = (NTimePoints, NVariables, NCells)".format(tube_id, datas.shape))
  # print(datas_dict)

  # optional slicing in time-dimension
  t_index_array = (times>=tplotmin) & (times<=tplotmax)
 
  # check if index array is empty
  if not np.any(t_index_array):
    raise IndexError('time index array is empty!')

  rho_data = datas[t_index_array, datas_dict['Density'], :]
  p_data = datas[t_index_array, datas_dict['Pressure'], :]
  
  if not 'PassiveScalars' in datas_dict:
   raise KeyError('No PassiveScalars data could be found!')
   
  rhoX_data = datas[t_index_array, datas_dict['PassiveScalars'], :]
  scalarX_data = rhoX_data / rho_data
  indexed_time = times[t_index_array]

  if not valueCheck(test_scalar, scalarX_data):
    raise ValueError("value Check failed!")
  

  print(scalarX_data.shape)
  print(np.min(scalarX_data), np.max(scalarX_data))

  time = []
  indices = []

  for step in range(scalarX_data.shape[0]):
    # print("Search for {}".format(test_scalar))
    current_time = indexed_time[step]
    if not valueCheck(test_scalar, scalarX_data[step,:]):
      continue
    
    idNeg, idPos = findNearest1D(scalarX_data[step,:], test_scalar)
    
    passive_scalar[0].append(interpolate1D(test_scalar, scalarX_data[step,:], p_data[step,:], (idNeg, idPos)))
    passive_scalar[1].append(interpolate1D(test_scalar, scalarX_data[step,:], rho_data[step,:], (idNeg, idPos)))
    
    time.append(current_time)
    indices.append(np.mean((idPos, idNeg)))
    
    if (idPos==scalarX_data.shape[1]-1) or (idNeg==scalarX_data.shape[1]-1):
      print("scalar has left the tube!")
      print("current time = {}".format(current_time))
      break

  passive_scalar = np.asarray_chkfinite(passive_scalar)

  time = np.array(time)
  pressure = passive_scalar[0]
  rho = passive_scalar[1]
  T = computeTemperature(pressure, rho)
  s = computeEntropie(T, rho, Dimless=Dimless)
  if not Dimless:
    rho *= ParameterNonDim['density']
    specVol = np.power(rho, -1)
    T *= ParameterNonDim['temperature']
  else:
    specVol = np.power(rho, -1)
    
  valueDict = {
    'time' : {
      'data' : time,
      'label': 'time',
      'symbol': 't',
      'dim': '[s]'
    },
    'id' : {
      'data' : indices,
      'label': 'id',
      'symbol': 'id',
      'dim': ''
    },
    'pressure' : {
      'data' : pressure,
      'label': 'pressure',
      'symbol': 'p',
      'dim': '[bar]'
    },
    'density' : {
      'data' : rho,
      'label': 'density',
      'symbol': 'rho',
      'dim': '[kg/m3]'
    },
    'temperature' : {
      'data' : T,
      'label': 'temperature',
      'symbol': 'T',
      'dim': '[K]'
    },
    'specificVolume' : {
      'data' : specVol,
      'label': 'specific Volume',
      'symbol': 'v',
      'dim': '[m3/kg]'
    },
    'specificEntropie' : {
      'data' : s,
      'label': 'specific Entropie',
      'symbol': 's',
      'dim': '[J/(kg*K)]'
    },
  }

  if Dimless:
    for key in valueDict.keys():
      if 'id' in key:
        continue
      valueDict[key]['dim'] = '[-]'

  # plot all variables over time
  for string in valueDict.keys():
    if 'time' in string:
      continue
    simplePlot('time', string, output_path)

  simplePlot('specificVolume', 'pressure', output_path)
  simplePlot('specificEntropie', 'temperature', output_path)
  simplePlot('time', 'id', output_path)

  f, ax = plt.subplots(nrows=1, ncols=1)
  x = specVol
  y = pressure
  ax.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], 
                scale_units='xy', angles='xy', scale=1)
  ax.plot(x[0], y[0], 'rx') #start
  ax.plot(x[-1], y[-1], 'gx') #end
  isoVol = np.linspace(specVol.min(), specVol.max(), 101 )
  const = 50.0 * 0.2 # p * v = const
  ax.plot(isoVol, const/isoVol, 'r--', label='isotherm')
  const = 50. * 0.2**gamma # p * v^gamma = const
  ax.plot(isoVol, const/(isoVol**gamma), 'b--', label='isentrop') 
  ax.set(xlabel="specific volume", ylabel="pressure")
  plt.legend(loc='best')
  f.savefig('{}/pvDiagrammArrow_tubeID{}.png'.format(output_path, tube_id), bbox_inches='tight', dpi=600)

  f, ax = plt.subplots(nrows=1, ncols=1)
  ax.quiver(s[:-1], T[:-1], s[1:]-s[:-1], T[1:]-T[:-1], 
                scale_units='xy', angles='xy', scale=1)
  ax.plot(s[0], T[0], 'rx')
  ax.plot(s[-1], T[-1], 'gx')
  ax.set(xlabel="specific entropie", ylabel="temperature")
  f.savefig('{}/tsDiagrammArrow_tubeID{}.png'.format(output_path, tube_id), bbox_inches='tight', dpi=600)

