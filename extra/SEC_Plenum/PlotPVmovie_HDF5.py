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
import matplotlib.gridspec as gridspec

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# check cli
if len(sys.argv)<3:
   errMsg = ('Not enough input arguments!\n'
               +'\t1. argument must be dataPath!\n'
               +'\t2. argument must be scalar number\n'
               +'\toptional argument is name of the inputfile\n'
               +'\toptional argument parallel working with module multiprocessing --parallel (default works on 8 CPUs)'
               +'\te.g. {} path 680 --config=inputfile.py --parallel=4'.format(sys.argv[0])
            )
   raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

DEFLAGRATION=False
if 'defl' in dataPath:
  DEFLAGRATION=True

# name of the inputfile is optional
optInputFilename = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--config=' in el ]
if not optInputFilename:
    optInputFilename = ['inputfile.py'] # default value 
inputfileName = optInputFilename[0]

PARALLEL=False
if any('--parallel' in arg for arg in sys.argv):
   from multiprocessing import Pool, cpu_count
   PARALLEL=True
   Num_CPUs = 8
   try:
      cpuString = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--parallel' in el ]
      if cpuString:
         Num_CPUs = cpuString[0]
   except: 
      pass
   if Num_CPUs>cpu_count():
      Num_CPUs = cpu_count()
   print(f'starting computations on {Num_CPUs} from {cpu_count()} available cores')

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
# plt.rcParams.update({'axes.axisbelow' : True})

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
test_scalar = float(sys.argv[2]) # e.g. 680.0
tplotmin = 200. #200.0
tplotmax = 400.0
Dimless = False # plot all variables dimless
# bool to read all existing HDF5 files

# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False

# what location / files should be used
INCLUDECONTROLSTATE = False
INCLUDETUBE = True
INCLUDETURBINEPLENUM = False

locDict = {
  'compressor' : 0, # control state
  'tube' : 1, # burning tube
  'turbinePlenum' : 2, # turbinePlenum
  'turbine' : 3, 
  'ambient' : 4
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
## list to collect all the data
## [0] --> pressure
## [1] --> density (later specific volume)
passive_scalar = [[],[]]
time = []
steps=[]
indices = []
location = []
tube_passive_scalar_times = {'start':None, 'stop':None}

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

  if not pvHelper.valueCheck(test_scalar, tube_scalarX):
    raise ValueError("value Check failed!")

  counter = 0
  FIRSTTIME = True
  for step in range(tube_scalarX.shape[0]):
    # print("Search for {}".format(test_scalar))
    current_time = indexed_time[step]
    if not pvHelper.valueCheck(test_scalar, tube_scalarX[step,:]):
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
      
    idNeg, idPos = pvHelper.findNearest1D(tube_scalarX[step,:], test_scalar)
    
    passive_scalar[0].append(pvHelper.interpolate1D(test_scalar, tube_scalarX[step,:], tube_p[step,:], (idNeg, idPos)))
    passive_scalar[1].append(pvHelper.interpolate1D(test_scalar, tube_scalarX[step,:], tube_rho[step,:], (idNeg, idPos)))
    time.append(current_time)
    steps.append(step)
    # indices.append(np.mean((idPos, idNeg)))
    indices.append(idPos)
    location.append(locDict['tube'])
    counter += 1

    if (idPos==tube_scalarX.shape[1]-1) or (idNeg==tube_scalarX.shape[1]-1):
      print("scalar has left the tube at current time = {}".format(current_time))
      tube_passive_scalar_times['stop'] = current_time
      break

if INCLUDETURBINEPLENUM:
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
    tube_scalarX = rhoX / rho

    if not pvHelper.valueCheck(test_scalar, tube_scalarX):
      continue
    
    pressure = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Pressure']], volume_fraction)
    
    # print(f'current time is = {current_time}')
    index = pvHelper.findNearest2D(tube_scalarX, test_scalar)

    passive_scalar[0].append(*pressure[index])
    passive_scalar[1].append(*rho[index])
    # find out why pressure[index] is a list??

    # print("current time = {}".format(current_time))
    # print(checkMax(test_scalar, tube_scalarX))
    time.append(current_time)
    location.append(locDict['turbinePlenum'])
    # indices.append(index) # index is 2d now!!

    if not pvHelper.checkMax(test_scalar, tube_scalarX):
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
valueDict['fuel']['data']=tube_fuelFraction

if Dimless:
  for key in valueDict.keys():
    if 'id' in key:
      continue
    valueDict[key]['dim'] = '[-]'

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9,2.5))
axs = axs.flatten()

print('plot quiver diagram')
pvHelper.quiverPlot('specificVolume', 'pressure', valueDict, axs[0], location, usetex=usetex)
pvHelper.quiverPlot('specificEntropy', 'temperature', valueDict, axs[1], location, usetex=usetex)

pvHelper.drawIsothermLine(50., 0.2, axs[0])
pvHelper.drawIsentropLine(50., 0.2, pv_class.gamma, axs[0])

# find the point when pressure is rising
p0, v0 = pvHelper.findPressureRisingPoint(valueDict['pressure']['data'], valueDict['specificVolume']['data'])
pvHelper.drawHugoniotCurve(p0, v0, pv_class.gamma, axs[0])
pvHelper.drawIsentropLine(p0, v0, pv_class.gamma, axs[0])

pvHelper.legend_without_duplicate_labels(axs[0])

ndig=3
if not Dimless:
  tube_passive_scalar_times = {key: value*pv_class.ref_time for key,value in tube_passive_scalar_times.items()}
  ndig = 5
fig.suptitle('passive scalar number {} was in tube [{}, {}]'.format(
                                  test_scalar,
                                  round(tube_passive_scalar_times['start'], ndig), 
                                  round(tube_passive_scalar_times['stop'], ndig) 
                                ))

figname = '{}/{}-quiver_{}-{}-Diagramm_tubeID{}'
if Dimless:
  figname += '_Dimless'
figname +='.png'
fig.savefig(figname.format(output_path, str(test_scalar).zfill(5), 'pv', 'Ts', tube_id), bbox_inches='tight', dpi=600)

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
        minT *= pv_class.ref_temperature
        maxT *= pv_class.ref_temperature
      
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
      if DEFLAGRATION:
        levels = np.linspace(1.0, 8., nlevels)
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
  newFolder = '{}/{}-timeseries'.format(output_path, test_scalar)
  os.makedirs(newFolder, exist_ok=True)

  plt.close('all')

  tube_extent, x, y_lower, y_upper = getTubeGeom(tubeHDF5)

  # plot in ax3 only upper half of the tube
  tube_extent_ax3 = tube_extent.copy()
  tube_extent_ax3[2] = 0.

  # plot in ax4 the lower part
  tube_extent_ax4 = tube_extent.copy()
  tube_extent_ax4[3] = 0.

  tube_temperature = pv_class.computeTemperature(tube_p[:,:], tube_rho[:,:], Dimless)
  
  print('Plotting timeseries')

  def plotFigure(i, step, valueDict, pv_class, PARALLEL=False):
    t = time[i]
    if PARALLEL:
      print(f'plotting picture {i}')
    
    sl = slice(i+2)

    fig = plt.figure(num=i, figsize=(9,10))
    
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

    pvHelper.quiverPlot('specificVolume', 'pressure', valueDict, ax1, location, sl=sl, usetex=usetex)
    pvHelper.quiverPlot('specificEntropy', 'temperature', valueDict, ax2, location, sl=sl, usetex=usetex)

    pvHelper.drawIsothermLine(50., 0.2, ax1)
    pvHelper.drawIsentropLine(50., 0.2, pv_class.gamma, ax1)
    
    p0, v0 = pvHelper.findPressureRisingPoint(valueDict['pressure']['data'], valueDict['specificVolume']['data'])
    pvHelper.drawHugoniotCurve(p0, v0, pv_class.gamma, ax1)
    pvHelper.drawIsentropLine(p0, v0, pv_class.gamma, ax1)
    
    pvHelper.legend_without_duplicate_labels(ax1)

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
    key = 'label'
    if usetex:
      key = 'latexLabel'
    cbarLabel = f"{valueDict['pressure'][key]} {valueDict['pressure']['dim']}"
    cbar.set_label(cbarLabel, rotation=270, labelpad=15)

    im_t = ax4.contourf(*temperature, extent=tube_extent_ax4, **props('temperature') )
    #ax4.fill_between(x, y_upper, np.max(y_upper), color='white')
    ax4.fill_between(x, y_lower, np.min(y_lower), color='white')
    ax4.axvline(x=x[indices[i]], color='k', ls='--')
    cbar = plt.colorbar(im_t, ax=ax4, shrink=0.9)
    cbarLabel = f"{valueDict['temperature'][key]} {valueDict['temperature']['dim']}"
    cbar.set_label(cbarLabel, rotation=270, labelpad=15)

    # print(f'fuel min = {np.min(fuel)}, max = {np.max(fuel)}')
    im_f = ax5.contourf(*fuel, extent=tube_extent, **props('fuel massfraction') )
    ax5.axvline(x=x[indices[i]], color='k', ls='--')
    cbar = plt.colorbar(im_f, ax=ax5, ticks=[0., 0.25, 0.5, 0.75, 1.0])
    cbarLabel = f"{valueDict['fuel'][key]} {valueDict['fuel']['dim']}"
    cbar.set_label(cbarLabel, rotation=270, labelpad=15)

    fig.suptitle('t = {:.8f} {}'.format(t, valueDict['time']['dim']), y=0.92, fontsize=14)
    fig.savefig('{}/{}-quiver_tubeID{}-{}.png'.format(newFolder, str(test_scalar).zfill(5), tube_id, str(i).zfill(5)), bbox_inches='tight', dpi=150)
    plt.close(i)

  if PARALLEL:
    values = ((i,step,valueDict,pv_class,PARALLEL) for i, step in enumerate(steps))

    with Pool(Num_CPUs) as pool:
      pool.starmap(plotFigure, values)
  else:

  # for i, (step, t) in enumerate(zip(steps, time)):
    for i, step in other.progressBar(steps, enumeration=True):
  # for i, step in enumerate(steps):
  #   ##print(f'step {i}')
      plotFigure(i, step, valueDict, pv_class)

  try:
    os.system('ffmpeg -framerate 10 -i {}/{}-quiver_tubeID{}-%5d.png -crf 20 {}/../{}Movie.mkv'.format(newFolder, str(test_scalar).zfill(5), tube_id, newFolder, str(test_scalar).zfill(5)) )
  except:
    print('ffmpeg could not be started')