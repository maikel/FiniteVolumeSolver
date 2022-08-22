import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
import amrex.other as other
import amrex.h5_data_processing as dataManip

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

other.import_file_as_module(os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref, L_ref, R_ref, R, gamma
from inputfile import D as diameter_tube

try:
  from inputfile import T_ref, n_tubes  
except:
  from inputfile import T_ref
  n_tubes=1

if not n_tubes==1:
  raise NotImplementedError('MultiTube Case not implemented yet!!')

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

def findNearest2D(array, value):
  absArray = np.abs(array-value)
  index = np.where( absArray == np.amin(absArray) )
  return index

def computeEntropy(T, rho, s0=0., Dimless=True):
  cv = getHeatCapacityConstantVolume(Dimless)
  Rspec = getSpecGasConstant(Dimless)
  s = np.zeros_like(rho)
  s[0] = s0 # set first element
  sdiff = (cv * np.log(T[1:]/T[:-1]) 
            + Rspec * np.log(rho[:-1]/rho[1:])
          )
  for i in range(s.shape[0]-1):
    s[i+1] = s[i]+sdiff[i]
  return s

def computeEnthalpy(T, Dimless=True):
  # H/m = h = cp * T
  cp = getHeatCapacityConstantPressure(Dimless)
  return cp*T

def computeTemperature(p, rho, Dimless=True):
  Rspec = getSpecGasConstant(Dimless)
  return p/(rho*Rspec)

def computeDensity(p, T, Dimless=True):
  Rspec = getSpecGasConstant(Dimless)
  return p/(T*Rspec)

def computePressure(rho, T, Dimless=True):
  Rspec = getSpecGasConstant(Dimless)
  return (rho*T*Rspec)

def checkMin(value, array):
  amin = np.amin(array)
  if (value < amin):
    print('scalar value: {} is lower than min scalar value: {}'.format(value, amin))
    return False
  else:
    return True

def checkMax(value, array):
  amax = np.amax(array)
  if (value > amax):
    print('scalar value: {} is greater than max scalar value: {}'.format(value, amax))
    return False
  else:
    return True


def valueCheck(value, array):
    MIN = checkMin(value, array)
    MAX = checkMax(value, array)
    if not MIN:
      return False
    if not MAX:
      return False
    return True

def interpolate1D(value, x, y, ids):
      data = np.interp( value, 
                [x[ids[0]], x[ids[1]]], 
                [y[ids[0]], y[ids[1]]] 
               )
      return data

def simplePlot(strX, strY, path, test_scalar, markerstyle='x', tube_id=0):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(valueDict[strX]['data'], valueDict[strY]['data'], 
              '-', marker=markerstyle)
    ax.set(
      xlabel=f"{valueDict[strX]['latexLabel']} {valueDict[strX]['dim']}", 
      ylabel=f"{valueDict[strY]['latexLabel']} {valueDict[strY]['dim']}"
      )
    fig.savefig('{}/{}-{}{}Diagramm_tubeID{}.png'.format(path, str(int(test_scalar)).zfill(5), valueDict[strY]['symbol'], valueDict[strX]['symbol'], tube_id), bbox_inches='tight')

def quiverPlot(strX, strY, path, test_scalar, location=[], tube_id=0):
  
  color = ['m', 'k', 'g']
  kwargs = {
    'scale_units' : 'xy', 
    'angles' : 'xy', 
    'scale' : 1
  }
  if location:
    kwargs['color'] = [color[i] for i in location]
  
  fig, ax = plt.subplots(nrows=1, ncols=1)
  x = valueDict[strX]['data']
  y = valueDict[strY]['data']
  
  ax.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], 
                **kwargs)
  # ax.plot(x[0], y[0], 'rx') #start
  # ax.plot(x[-1], y[-1], 'gx') #end

  if 'pressure' in valueDict[strY]['label']:
    isoVol = np.linspace(x.min(), x.max(), 101 )
    const = 50.0 * 0.2 # p * v = const
    ax.plot(isoVol, const/isoVol, 'r--', label='isotherm')
    const = 50. * 0.2**gamma # p * v^gamma = const
    ax.plot(isoVol, const/(isoVol**gamma), 'b--', label='isentrop') 
    plt.legend(loc='best')

  ax.set(
      xlabel=f"{valueDict[strX]['latexLabel']} {valueDict[strX]['dim']}", 
      ylabel=f"{valueDict[strY]['latexLabel']} {valueDict[strY]['dim']}"
      )
  fig.savefig('{}/{}-quiver_{}{}Diagramm_tubeID{}.png'.format(path, str(int(test_scalar)).zfill(5), valueDict[strY]['symbol'], valueDict[strX]['symbol'], tube_id), bbox_inches='tight', dpi=600)

def getControlStateData(timepoint, csData, csTimes):
  return csData[csTimes>=timepoint, :][0]

#-------------------------------------------------

plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization/PVDiagramm/'.format(outPath)

os.makedirs(output_path, exist_ok=True)

#---------------------------
test_scalar = 680.0 # initial passive scalar value
tplotmin = 300. #200.0
tplotmax = 400.0
Dimless = True
# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False
INCLUDECONTROLSTATE = True
INCLUDETUBE = True
INCLUDEPLENUM = True

locDict = {
  'compressor' : 0, # control state
  'tube' : 1, # burning tube
  'turbine' : 2 # turbine plenum
}
#-------------------------------------------
# load ControlState data

if INCLUDECONTROLSTATE:
  filename_basic = "{}/ControlState.h5".format(dataPath)
  print("Read in data from {}".format(filename_basic))

  if RESTARTEDSIMULATION:
    csData, csTimes, csDict = io.h5_load_restartedTimeseries(filename_basic)
  else:
    csData, csTimes, csDict = io.h5_load_timeseries(filename_basic)

#-------------------------------------------
# ambientPressure = 1. # 1 bar
# ambientTemperature = 273.15 / T_ref # 0°C
# ambientDensity = computeDensity(ambientPressure, ambientTemperature)
# passive_scalar = [[ambientPressure],[ambientDensity]]
# time = [tplotmin]


## list to collect all the data
## [0] --> pressure
## [1] --> density (later specific volume)
passive_scalar = [[],[]]
time = []
indices = []
location = []

if INCLUDETUBE:
  tube_id = 0
  filename_basic = '{}/Tube{}.h5'.format(dataPath, tube_id)
  # datas, times, datas_dict = io.h5_load_timeseries(filename_basic)
  extent_1d = io.h5_load_get_extent_1D(filename_basic)

  if RESTARTEDSIMULATION:
    datas, times, datas_dict = io.h5_load_restartedTimeseries(filename_basic)
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

  counter = 0
  for step in range(scalarX_data.shape[0]):
    # print("Search for {}".format(test_scalar))
    current_time = indexed_time[step]
    if not valueCheck(test_scalar, scalarX_data[step,:]):
      continue
    
    if (counter==0) and INCLUDECONTROLSTATE:
      csdata = getControlStateData(current_time, csData, csTimes)
      # 'compressor_pressure', compressor_temperature'
      passive_scalar[0].append(csdata[csDict['compressor_pressure']])
      csDensity = computeDensity(
            csdata[csDict['compressor_pressure']],
            csdata[csDict['compressor_temperature']] )
      passive_scalar[1].append(csDensity)
      time.append(current_time)
      location.append(locDict['compressor'])
      del csDensity, csdata
      
    idNeg, idPos = findNearest1D(scalarX_data[step,:], test_scalar)
    
    passive_scalar[0].append(interpolate1D(test_scalar, scalarX_data[step,:], p_data[step,:], (idNeg, idPos)))
    passive_scalar[1].append(interpolate1D(test_scalar, scalarX_data[step,:], rho_data[step,:], (idNeg, idPos)))
    time.append(current_time)
    indices.append(np.mean((idPos, idNeg)))
    location.append(locDict['tube'])
    counter += 1

    if (idPos==scalarX_data.shape[1]-1) or (idNeg==scalarX_data.shape[1]-1):
      print("scalar has left the tube!")
      print("current time = {}".format(current_time))
      break

if INCLUDEPLENUM:
  # Get times and nSteps from plenum
  plenum_times = io.h5_load_timepoints(plenum)
  nSteps = plenum_times.shape[0]
  t_bool_array = (plenum_times>=time[-1]) & (plenum_times<=tplotmax)
  t_index_array = np.nonzero(t_bool_array)[0] # returns tuple: (array([10 11 ...]), )

  # check if index array is empty
  if not np.any(t_index_array):
    raise IndexError('time index array is empty!')

  # for i in other.progressBar(t_index_array):
  for i in t_index_array:
    plenum_variables = ["Pressure", "Density", "PassiveScalars", 'vfrac']
    
    plenum_data, current_time, plenum_extent, plenum_dict = io.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
    
    volume_fraction = plenum_data[plenum_dict['vfrac']]
    rhoX = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['PassiveScalars']], volume_fraction)
    rho = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Density']], volume_fraction)
    scalarX_data = rhoX / rho

    if not valueCheck(test_scalar, scalarX_data):
      continue
    
    pressure = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Pressure']], volume_fraction)
    
    # print(f'current time is = {current_time}')
    index = findNearest2D(scalarX_data, test_scalar)

    # possible interpolation...
    passive_scalar[0].append(*pressure[index])
    passive_scalar[1].append(*rho[index])
    # find out why pressure[index] is a list??

    # print("current time = {}".format(current_time))
    # print(checkMax(test_scalar, scalarX_data))
    time.append(current_time)
    location.append(locDict['turbine'])
    # indices.append(index) # index is 2d now!!

    if not checkMax(test_scalar, scalarX_data):
      print("scalar has left the plenum!")
      print("current time = {}".format(current_time))
      break
    del pressure, rho


##-------------------------------------------------------
## arange data

## add ambient data
# passive_scalar[0].append(ambientPressure)
# passive_scalar[1].append(ambientDensity)
# time.append(current_time)

passive_scalar = np.asarray_chkfinite(passive_scalar)

time = np.array(time)
pressure = passive_scalar[0]
rho = passive_scalar[1]
T = computeTemperature(pressure, rho)
s = computeEntropy(T, rho, Dimless=Dimless)
h = computeEnthalpy(T, Dimless)
if not Dimless:
  rho *= ParameterNonDim['density']
  specVol = np.power(rho, -1)
  T *= ParameterNonDim['temperature']
  time *= ParameterNonDim['time']
else:
  specVol = np.power(rho, -1)
  
valueDict = {
  'time' : {
    'data' : time,
    'label': 'time',
    'symbol': 't',
    'latexLabel': r'$t$',
    'dim': '[s]'
  },
  # 'id' : {
  #   'data' : indices,
  #   'label': 'id',
  #   'symbol': 'id',
  #   'dim': ''
  # },
  'pressure' : {
    'data' : pressure,
    'label': 'pressure',
    'symbol': 'p',
    'latexLabel': r'$p$',
    'dim': '[bar]'
  },
  'density' : {
    'data' : rho,
    'label': 'density',
    'symbol': 'rho',
    'latexLabel': r'$\rho$',
    'dim': '[kg/m3]'
  },
  'temperature' : {
    'data' : T,
    'label': 'temperature',
    'symbol': 'T',
    'latexLabel': r'$T$',
    'dim': '[K]'
  },
  'specificVolume' : {
    'data' : specVol,
    'label': 'specific Volume',
    'symbol': 'v',
    'latexLabel': r'$v$',
    'dim': '[m3/kg]'
  },
  'specificEntropy' : {
    'data' : s,
    'label': 'specific Entropy difference',
    'symbol': 's',
    'latexLabel': r'$\Delta s$',
    'dim': '[J/(kg*K)]'
  },
  'specificEnthalpy' : {
    'data' : h,
    'label': 'specific Enthalpy',
    'symbol': 'h',
    'latexLabel': r'$h$',
    'dim': '[J/kg]'
  },
}

if Dimless:
  for key in valueDict.keys():
    if 'id' in key:
      continue
    valueDict[key]['dim'] = '[-]'

# # plot all variables over time
# for string in valueDict.keys():
#   if 'time' in string:
#     continue
#   simplePlot('time', string, output_path, test_scalar)

# simplePlot('specificVolume', 'pressure', output_path, test_scalar)
# simplePlot('specificEntropy', 'temperature', output_path, test_scalar)
# simplePlot('specificEntropy', 'specificEnthalpy', output_path, test_scalar)
  
quiverPlot('specificVolume', 'pressure', output_path, test_scalar, location)
quiverPlot('specificEntropy', 'temperature', output_path, test_scalar, location)
quiverPlot('specificEntropy', 'specificEnthalpy', output_path, test_scalar, location)