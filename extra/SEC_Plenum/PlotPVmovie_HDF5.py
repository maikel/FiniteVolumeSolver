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
import matplotlib.gridspec as gridspec

import h5py
import math

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
# plt.rcParams.update({'axes.axisbelow' : True})

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
  Rspec = getSpecGasConstant(True)
  if Dimless:
    return p/(rho*Rspec)
  else:
    return p/(rho*Rspec)*ParameterNonDim['temperature']

def computeDensity(p, T, Dimless=True):
  Rspec = getSpecGasConstant(True)
  if Dimless:
    return p/(T*Rspec)
  else: 
    return p/(T*Rspec)*ParameterNonDim['density']

def computePressure(rho, T, Dimless=True):
  Rspec = getSpecGasConstant(True)
  if Dimless:
    return (rho*T*Rspec)
  else:
    return (rho*T*Rspec)*ParameterNonDim['pressure']

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

def quiverPlot(strX, strY, ax, location=[], tube_id=0, sl=()):
  print(f'Plotting {strY}-{strX} diagram')
  color = ['m', 'k', 'g']
  kwargs = {
    'scale_units' : 'xy', 
    'angles' : 'xy', 
    'scale' : 1
  }
  if location:
    kwargs['color'] = [color[i] for i in location]
  
  x = valueDict[strX]['data']
  y = valueDict[strY]['data']

  ax.set(
    xlim = (math.floor(1.1*x.min()), round(1.1*x.max(), 1)),
    ylim = (math.floor(1.1*y.min()), round(1.1*y.max()))
  )
  
  # ax.plot(x[0], y[0], 'rx') #start
  # ax.plot(x[-1], y[-1], 'gx') #end

  if 'pressure' in valueDict[strY]['label']:
    isoVol = np.linspace(x.min(), x.max(), 101 )
    const = 50.0 * 0.2 # p * v = const
    ax.plot(isoVol, const/isoVol, 'r--', label='isotherm '+r'$p\cdot v$')
    const = 50. * 0.2**gamma # p * v^gamma = const
    ax.plot(isoVol, const/(isoVol**gamma), 'b--', label='isentrop '+r'$p\cdot v^\gamma$') 
    ax.legend(loc='best')
    const = 10. * 0.2**gamma # p * v^gamma = const
    ax.plot(isoVol, const/(isoVol**gamma), 'b--') 
    

  if sl:
    x = x[sl]
    y = y[sl]
  
  ax.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], 
                **kwargs)

  ax.set(
      xlabel=f"{valueDict[strX]['latexLabel']} {valueDict[strX]['dim']}", 
      ylabel=f"{valueDict[strY]['latexLabel']} {valueDict[strY]['dim']}"
      )

def getControlStateData(timepoint, csData, csTimes):
  return csData[csTimes>=timepoint, :][0]

#-------------------------------------------------

plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization/PVDiagramm/'.format(outPath)

os.makedirs(output_path, exist_ok=True)

#---------------------------
# initial passive scalar value
try:
  test_scalar = float(sys.argv[3])
except: 
  test_scalar = 680.0
tplotmin = 300. #200.0
tplotmax = 400.0
Dimless = False # plot all variables dimless
# bool to read all existing HDF5 files

# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False

# what location / files should be used
INCLUDECONTROLSTATE = False
INCLUDETUBE = True
INCLUDEPLENUM = False

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
steps=[]
indices = []
location = []

if INCLUDETUBE:
  tube_id = 0
  tubeHDF5 = '{}/Tube{}.h5'.format(dataPath, tube_id)
  # datas, times, datas_dict = io.h5_load_timeseries(tubeHDF5)
  extent_1d = io.h5_load_get_extent_1D(tubeHDF5)

  if RESTARTEDSIMULATION:
    datas, times, datas_dict = io.h5_load_restartedTimeseries(tubeHDF5)
  else:
    print("[Tube{}] Read in data from {}".format(tube_id, tubeHDF5))
    datas, times, datas_dict = io.h5_load_timeseries(tubeHDF5)


  datas = np.squeeze(datas) # remove last axis
  print("[Tube{}] data shape from tube is {} = (NTimePoints, NVariables, NCells)".format(tube_id, datas.shape))
  # print(datas_dict)

  # optional slicing in time-dimension
  t_index_array = (times>=tplotmin) & (times<=tplotmax)

  # check if index array is empty
  if not np.any(t_index_array):
    raise IndexError('time index array is empty!')

  tube_rho = datas[t_index_array, datas_dict['Density'], :]
  tube_p = datas[t_index_array, datas_dict['Pressure'], :]
  tube_fuelFraction = datas[t_index_array, datas_dict['Species'], :] / tube_rho

  if not 'PassiveScalars' in datas_dict:
    raise KeyError('No PassiveScalars data could be found!')
    
  rhoX_data = datas[t_index_array, datas_dict['PassiveScalars'], :]
  tube_scalarX = rhoX_data / tube_rho
  indexed_time = times[t_index_array]

  if not valueCheck(test_scalar, tube_scalarX):
    raise ValueError("value Check failed!")


  print(tube_scalarX.shape)
  print(np.min(tube_scalarX), np.max(tube_scalarX))

  counter = 0
  for step in range(tube_scalarX.shape[0]):
    # print("Search for {}".format(test_scalar))
    current_time = indexed_time[step]
    if not valueCheck(test_scalar, tube_scalarX[step,:]):
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
      
    idNeg, idPos = findNearest1D(tube_scalarX[step,:], test_scalar)
    
    passive_scalar[0].append(interpolate1D(test_scalar, tube_scalarX[step,:], tube_p[step,:], (idNeg, idPos)))
    passive_scalar[1].append(interpolate1D(test_scalar, tube_scalarX[step,:], tube_rho[step,:], (idNeg, idPos)))
    time.append(current_time)
    steps.append(step)
    # indices.append(np.mean((idPos, idNeg)))
    indices.append(idPos)
    location.append(locDict['tube'])
    counter += 1

    if (idPos==tube_scalarX.shape[1]-1) or (idNeg==tube_scalarX.shape[1]-1):
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
    tube_scalarX = rhoX / rho

    if not valueCheck(test_scalar, tube_scalarX):
      continue
    
    pressure = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Pressure']], volume_fraction)
    
    # print(f'current time is = {current_time}')
    index = findNearest2D(tube_scalarX, test_scalar)

    # possible interpolation...
    passive_scalar[0].append(*pressure[index])
    passive_scalar[1].append(*rho[index])
    # find out why pressure[index] is a list??

    # print("current time = {}".format(current_time))
    # print(checkMax(test_scalar, tube_scalarX))
    time.append(current_time)
    location.append(locDict['turbine'])
    # indices.append(index) # index is 2d now!!

    if not checkMax(test_scalar, tube_scalarX):
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
T = computeTemperature(pressure, rho, Dimless)
s = computeEntropy(T, rho, Dimless=Dimless)
h = computeEnthalpy(T, Dimless)
if not Dimless:
  rho *= ParameterNonDim['density']
  specVol = np.power(rho, -1)
  #T *= ParameterNonDim['temperature']
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
    'dim': r'$[kg/m^3]$'
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
    'dim': r'$[m^3/kg]$'
  },
  'specificEntropy' : {
    'data' : s,
    'label': 'specific Entropy difference',
    'symbol': 's',
    'latexLabel': r'$\Delta s$',
    'dim': r'$[J/(kg\cdot K)]$'
  },
  'specificEnthalpy' : {
    'data' : h,
    'label': 'specific Enthalpy',
    'symbol': 'h',
    'latexLabel': r'$h$',
    'dim': '[J/kg]'
  },
  'fuel' : {
    'data' : tube_fuelFraction,
    'label': 'fuel',
    'symbol': 'Yf',
    'latexLabel': r'$Y_F$',
    'dim': '[-]'
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

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9,2.5))
axs = axs.flatten()

quiverPlot('specificVolume', 'pressure', axs[0], location)
quiverPlot('specificEntropy', 'temperature', axs[1], location)

fig.savefig('{}/{}-quiver_{}-{}-Diagramm_tubeID{}.png'.format(output_path, str(int(test_scalar)).zfill(5), 'pv', 'Ts', tube_id), bbox_inches='tight', dpi=600)


# quiverPlot('specificVolume', 'pressure', output_path, test_scalar, location)
# quiverPlot('specificEntropy', 'temperature', output_path, test_scalar, location)
# quiverPlot('specificEntropy', 'specificEnthalpy', output_path, test_scalar, location)

#---------------------------------------------------------------
def getTubeGeom(filename):
  tube_variables = ["Pressure"] #, "PassiveScalars", "Density"]
  _, _, tube_extent, _ = io.h5_load_spec_timepoint_variable(filename, 0, tube_variables)

  midpoint = 0.5 * (tube_extent[2] + tube_extent[3])
  x = np.linspace(tube_extent[0], tube_extent[1], tube_n_cells)
  y_upper = midpoint + 0.5 * np.array([Area(xi) for xi in x]) * diameter_tube
  y_lower = midpoint - 0.5 * np.array([Area(xi) for xi in x]) * diameter_tube
  # y_upper_scale = 0.1
  # y_upper += y_upper_scale*np.max(y_upper) # increase y_upper by y_upper_scale (in percent) to avoid overlapping in the image
  lower = midpoint - 2.0 * diameter_tube
  upper = midpoint + 2.0 * diameter_tube
  tube_extent[2] = lower
  tube_extent[3] = upper
  
  return tube_extent, x, y_lower, y_upper

def stackTubeDataTo2D(Tube_datalist):
    # all Tubedata is 1D but for contourf we need at least 2D data. so simply stack twice the 1d array
    for i, el in enumerate(Tube_datalist):
      el = np.squeeze(el)
      Tube_datalist[i] = np.stack((el,el))
    return Tube_datalist

def props(title):
   nlevels = 96
   props = {
      'extend': 'max',
   }
   if  title == 'temperature':
      minT = 1.
      maxT = 8. #11.
      if not Dimless:
        minT *= ParameterNonDim['temperature']
        maxT *= ParameterNonDim['temperature']
      
      levels = np.linspace(minT, maxT, nlevels )
      props.update(
         {
            'levels': levels,
            'vmin': levels[0],
            'vmax': levels[-1],
            'cmap': 'coolwarm' # 'afmhot_r'
         }
      )
   elif title == 'pressure':
      levels = np.linspace(1.0, 20., nlevels)
      props.update(
         {
         'levels': levels,
         'vmin': levels[0],
         'vmax': levels[-1],
         'cmap': 'jet'
         }
      )
   elif title == 'fuel massfraction':
      levels = np.linspace(0., 1., nlevels)
      props = {
         'levels': levels,
         'vmin': levels[0],
         'vmax': levels[-1],
         'cmap': 'binary'
         }
   return props
#----------------------------------------------------

if INCLUDETUBE:
  newFolder = '{}/{}-timeseries'.format(output_path, int(test_scalar))
  os.makedirs(newFolder, exist_ok=True)

  plt.close('all')

  tube_extent, x, y_lower, y_upper = getTubeGeom(tubeHDF5)

  # plot in ax3 only upper half of the tube
  tube_extent_ax3 = tube_extent.copy()
  tube_extent_ax3[2] = 0.

  # plot in ax4 the lower part
  tube_extent_ax4 = tube_extent.copy()
  tube_extent_ax4[3] = 0.

  tube_temperature = computeTemperature(tube_p[:,:], tube_rho[:,:], Dimless)

  for i, (step, t) in enumerate(zip(steps, time)):
    print(f'step {i}')
    sl = slice(i+2)

    fig = plt.figure(figsize=(9,10))
    
    # create gridspec in figure
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[6,6,1])
    
    # create nested gridspec for two quiverplots
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0])

    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[0, 1])

    # create another nested gridspec for combustion tube plot
    gs01 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[1], hspace=0.)

    ax3 = fig.add_subplot(gs01[0, :])
    ax4 = fig.add_subplot(gs01[1, :])

    ax5 = fig.add_subplot(gs0[2])

    quiverPlot('specificVolume', 'pressure', ax1, location, sl=sl)
    quiverPlot('specificEntropy', 'temperature', ax2, location, sl=sl)

    #print(f'pressure min={pressure.min()} and max={pressure.max()}')
    pressure = stackTubeDataTo2D([tube_p[step,:]])
    temperature = stackTubeDataTo2D([ tube_temperature[step,:] ])
    fuel = stackTubeDataTo2D([ tube_fuelFraction[step,:] ])

    im_p = ax3.contourf(*pressure, extent=tube_extent_ax3, **props('pressure') )
    ax3.fill_between(x, y_upper, np.max(y_upper), color='white')
    #ax3.fill_between(x, y_lower, np.min(y_lower), color='white')
    ax3.axvline(x=x[indices[i]], color='k', ls='--')
    ax3.xaxis.set_ticklabels([]) #hide only the labels
    cbar = plt.colorbar(im_p, ax=ax3, shrink=0.9)
    cbarLabel = f"{valueDict['pressure']['latexLabel']} {valueDict['pressure']['dim']}"
    cbar.set_label(cbarLabel, rotation=270, labelpad=15)

    im_t = ax4.contourf(*temperature, extent=tube_extent_ax4, **props('temperature') )
    #ax4.fill_between(x, y_upper, np.max(y_upper), color='white')
    ax4.fill_between(x, y_lower, np.min(y_lower), color='white')
    ax4.axvline(x=x[indices[i]], color='k', ls='--')
    cbar = plt.colorbar(im_t, ax=ax4, shrink=0.9)
    cbarLabel = f"{valueDict['temperature']['latexLabel']} {valueDict['temperature']['dim']}"
    cbar.set_label(cbarLabel, rotation=270, labelpad=15)

    # print(f'fuel min = {np.min(fuel)}, max = {np.max(fuel)}')
    im_f = ax5.contourf(*fuel, extent=tube_extent, **props('fuel massfraction') )
    ax5.axvline(x=x[indices[i]], color='k', ls='--')
    cbar = plt.colorbar(im_f, ax=ax5, ticks=[0., 0.25, 0.5, 0.75, 1.0])
    cbarLabel = f"{valueDict['fuel']['latexLabel']} {valueDict['fuel']['dim']}"
    cbar.set_label(cbarLabel, rotation=270, labelpad=15)

    fig.suptitle('t = {:.8f} {}'.format(t, valueDict['time']['dim']), y=0.92, fontsize=14)
    fig.savefig('{}/{}-quiver_tubeID{}-{}.png'.format(newFolder, str(int(test_scalar)).zfill(5), tube_id, str(i).zfill(5)), bbox_inches='tight', dpi=150)
    plt.close('all')

  os.system('ffmpeg -framerate 10 -i {}/{}-quiver_tubeID{}-%5d.png -crf 20 {}/../{}Movie.mkv'.format(newFolder, str(int(test_scalar)).zfill(5), tube_id, newFolder, str(int(test_scalar)).zfill(5)) )