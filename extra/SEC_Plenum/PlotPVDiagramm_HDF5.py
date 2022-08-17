import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
import amrex.other as other
import amrex.h5_data_processing as dataManip
import amrex.pvDiagrammHelper as pvHelper

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import json

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# check cli
if len(sys.argv)<3:
   errMsg = ('Not enough input arguments!\n'
               +'\t1. argument must be dataPath!\n'
               +'\t2. argument must be list of scalar numbers without white spaces\n'
               +'\toptional argument is name of the inputfile\n'
               +'\te.g. {} path [680, 690, ...] --config=inputfile.py'.format(sys.argv[0])
            )
   raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

# name of the inputfile is optional
optional = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--config=' in el ]
if not optional:
    optional = ['inputfile.py'] # default value 
inputfileName = optional[0]

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

usetex = matplotlib.checkdep_usetex(True)
if usetex:
  # Use LaTeX to write all text
  tex_fonts.update({"text.usetex": True, 
                  "font.family": "serif"}) 
plt.rcParams.update(tex_fonts)
plt.rcParams.update({'axes.grid' : True})
plt.rcParams.update({'axes.axisbelow' : True})

#--------------------------------------------------------------------------
print('initiate Attribute class')
pv_class = pvHelper.pv_diagramm_class()
pv_class.setAttribute('specGasConstant', R) # [-]
pv_class.setAttribute('gamma', gamma) # [-]
pv_class.setAttribute('ref_time', t_ref) # [s]
pv_class.setAttribute('ref_length', L_ref) # [m]
pv_class.setAttribute('ref_pressure', p_ref) # [Pa]
pv_class.setAttribute('ref_density' , rho_ref) # [kg/m^3]
pv_class.setAttribute('ref_temperature', T_ref) # [K]
pv_class.setAttribute('ref_velocity', u_ref) # [m/s]
pv_class.setAttribute('ref_GasConstant', R_ref) # [J/(kg*K)]

print('Attribute class contains following members:')
print(pv_class.prettyprint())

#-------------------------------------------------

outPath = dataPath
output_path = '{}/Visualization/PVDiagramm/'.format(outPath)

os.makedirs(output_path, exist_ok=True)

#---------------------------
# initial passive scalar value

temp = sys.argv[2].strip('[]()')
test_scalar_list = [pvHelper.convertStringToFloat(el) for el in temp.split(',')] # e.g.[680,690]
print('Following passive scalar numbers could be parsed: {}'.format(test_scalar_list))

tplotmin = 200. #200.0
tplotmax = 400.0
Dimless = True
# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False
INCLUDECONTROLSTATE = False
INCLUDETUBE = True
INCLUDEPLENUM = False

locDict = {
  'compressor' : 0, # control state
  'tube' : 1, # burning tube
  'turbine' : 2 # turbine plenum
}
#-------------------------------------------

for test_scalar in test_scalar_list:
  print('\nparsing data for passive scalar number {}'.format(test_scalar))
  ## list to collect all the data
  ## [0] --> pressure
  ## [1] --> density (later specific volume)
  passive_scalar = [[],[]]
  time = []
  indices = []
  location = []
  tube_passive_scalar_times = {'start':None, 'stop':None}
  
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
  # load Tube0 data
  if INCLUDETUBE:
    tube_id = 0
    filename_basic = '{}/Tube{}.h5'.format(dataPath, tube_id)
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

    if not pvHelper.valueCheck(test_scalar, scalarX_data):
      raise ValueError("value Check failed!")

    # print(scalarX_data.shape)
    # print(np.min(scalarX_data), np.max(scalarX_data))

    counter = 0
    FIRSTTIME = True
    for step in range(scalarX_data.shape[0]):
      # print("Search for {}".format(test_scalar))
      current_time = indexed_time[step]
      if not pvHelper.valueCheck(test_scalar, scalarX_data[step,:]):
        continue
      else:
        if FIRSTTIME:
          print("scalar has entered the tube at current time = {}".format(current_time))
          tube_passive_scalar_times['start'] = current_time
          FIRSTTIME = False

      if (counter==0) and INCLUDECONTROLSTATE:
        csdata = pvHelper.getControlStateData(current_time, csData, csTimes)
        # 'compressor_pressure', compressor_temperature'
        passive_scalar[0].append(csdata[csDict['compressor_pressure']])
        csDensity = pv_class.computeDensity(
              csdata[csDict['compressor_pressure']],
              csdata[csDict['compressor_temperature']] )
        passive_scalar[1].append(csDensity)
        time.append(current_time)
        location.append(locDict['compressor'])
        del csDensity, csdata
      
      idNeg, idPos = pvHelper.findNearest1D(scalarX_data[step,:], test_scalar)
      
      passive_scalar[0].append(pvHelper.interpolate1D(test_scalar, scalarX_data[step,:], p_data[step,:], (idNeg, idPos)))
      passive_scalar[1].append(pvHelper.interpolate1D(test_scalar, scalarX_data[step,:], rho_data[step,:], (idNeg, idPos)))
      time.append(current_time)
      indices.append(np.mean((idPos, idNeg)))
      location.append(locDict['tube'])
      counter += 1

      if (idPos==scalarX_data.shape[1]-1) or (idNeg==scalarX_data.shape[1]-1):
        print("scalar has left the tube at current time = {}".format(current_time))
        tube_passive_scalar_times['stop'] = current_time
        break

  #------------------------------
  # load Plenum data
  if INCLUDEPLENUM:
    plenum = "{}/Plenum.h5".format(dataPath)
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

      if not pvHelper.valueCheck(test_scalar, scalarX_data):
        continue
      
      pressure = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Pressure']], volume_fraction)
      
      # print(f'current time is = {current_time}')
      index = pvHelper.findNearest2D(scalarX_data, test_scalar)

      # possible interpolation...
      passive_scalar[0].append(*pressure[index])
      passive_scalar[1].append(*rho[index])
      # find out why pressure[index] is a list??

      # print("current time = {}".format(current_time))
      # print(checkMax(test_scalar, scalarX_data))
      time.append(current_time)
      location.append(locDict['turbine'])
      # indices.append(index) # index is 2d now!!

      if not pvHelper.checkMax(test_scalar, scalarX_data):
        print("scalar has left the plenum!")
        print("current time = {}".format(current_time))
        break
      del pressure, rho


  ##-------------------------------------------------------
  ## arange data
  passive_scalar = np.asarray_chkfinite(passive_scalar)

  time = np.array(time)
  pressure = passive_scalar[0]
  rho = passive_scalar[1]
  T = pv_class.computeTemperature(pressure, rho, Dimless)
  s = pv_class.computeEntropy(T, rho, Dimless)
  h = pv_class.computeEnthalpy(T, Dimless)
  if not Dimless:
    rho *= pv_class.ref_density
    specVol = np.power(rho, -1)
    time *= pv_class.ref_time
  else:
    specVol = np.power(rho, -1)
    
  from amrex.pvDiagrammHelper import valueDict
  valueDict['time']['data']=time
  valueDict['pressure']['data']=pressure
  valueDict['density']['data']=rho
  valueDict['temperature']['data']=T
  valueDict['specificVolume']['data']=specVol
  valueDict['specificEntropy']['data']=s
  valueDict['specificEnthalpy']['data']=h

  if Dimless:
    for key in valueDict.keys():
      if 'id' in key:
        continue
      valueDict[key]['dim'] = '[-]'


  # simplePlot('specificVolume', 'pressure', output_path, test_scalar)
  # simplePlot('specificEntropy', 'temperature', output_path, test_scalar)
  # simplePlot('specificEntropy', 'specificEnthalpy', output_path, test_scalar)
  
  fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9,2.5))
  axs = axs.flatten()

  print('plot quiver diagram')
  pvHelper.quiverPlot('specificVolume', 'pressure', valueDict, axs[0], location, usetex=usetex)
  pvHelper.quiverPlot('specificEntropy', 'temperature', valueDict, axs[1], location, usetex=usetex)

  pvHelper.drawIsothermLine(50., 0.2, axs[0])
  pvHelper.drawIsentropLine(50., 0.2, pv_class.gamma, axs[0])
  
  p0, v0 = pvHelper.findPressureRisingPoint(valueDict['pressure']['data'], valueDict['specificVolume']['data'])
  pvHelper.drawHugoniotCurve(p0, v0, pv_class.gamma, axs[0])
  pvHelper.drawIsentropLine(p0, v0, pv_class.gamma, axs[0])
  
  pvHelper.legend_without_duplicate_labels(axs[0])

  fig.suptitle('passive scalar number {} was in tube [{}, {}]'.format(
                                  test_scalar,
                                  round(tube_passive_scalar_times['start'], 2), 
                                  round(tube_passive_scalar_times['stop'], 2) 
                                ))
  
  figname = '{}/{}-quiver_{}-{}-Diagramm_tubeID{}'
  if Dimless:
    figname += '_Dimless'
  figname +='.png'
  fig.savefig(figname.format(output_path, str(int(test_scalar)).zfill(5), 'pv', 'Ts', tube_id), bbox_inches='tight', dpi=600)