import sys, os
# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
import amrex.h5_data_processing as dataManip
from amrex.other import import_file_as_module

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

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

output_path = '{}/Visualization'.format(dataPath)

# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False

try:
  inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
  inputfileName = 'SEC_Plenum_Arrhenius.py'

import_file_as_module(os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import t_ref, T_ref, ControlOptions

try:
  from inputfile import T_ref, n_tubes
except:
  from inputfile import T_ref
  n_tubes=1
  

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

#plt.style.use('seaborn')
tex_fonts = {
    # Use LaTeX to write all text
    # "text.usetex": True,
    # "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 12,# 10,
    "axes.titlesize": 12,# 10,
    "axes.labelsize": 12,# 10,
    "font.size": 12,# 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 12, # 10,
    "xtick.labelsize": 10,# 8,
    "ytick.labelsize": 10,# 8
}

CROPPED=True

plt.rcParams.update(tex_fonts)
plt.rcParams.update({'axes.grid' : False})

os.makedirs(output_path, exist_ok=True)


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
  tplotmin = 398 #395.0
  tplotmax = 400.0
  del datas_dict['PassiveScalars']

  t_index_array = (times>=tplotmin) & (times<=tplotmax)
  # check if index array is empty
  if not np.any(t_index_array):
    raise IndexError('time index array is empty!')

  # t_index_array = (times>=tplotmin)
  
  rho_data = datas[t_index_array, datas_dict['Density'], :]
  rhou_data = datas[t_index_array, datas_dict['Momentum'], :]
  rhoE_data = datas[t_index_array, datas_dict['Energy'], :]
  rhoF_data = datas[t_index_array, datas_dict['Species'], :]
  p_data = datas[t_index_array, datas_dict['Pressure'], :]
  c_data = datas[t_index_array, datas_dict['SpeedOfSound'], :]
  if 'PassiveScalars' in datas_dict:
    rhoX_data = datas[t_index_array, datas_dict['PassiveScalars'], :]
    X = rhoX_data / rho_data
  T_data = p_data / rho_data
  F_data = rhoF_data / rho_data
  Ma = rhou_data / rho_data / c_data

  x0 = extent_1d[0]
  xEnd = extent_1d[1]
  t0 = times[t_index_array][0]
  tEnd = round( times[t_index_array][-1], 2) # round tEnd to display the last yticklabel
  print("[Tube{}] tEnd is {}".format(tube_id, tEnd))

  # print out the first occurence of min/max value 
  dataManip.printSimpleStatsTubeData(p_data, 'Pressure', times[t_index_array], tube_id)
  dataManip.printSimpleStatsTubeData(T_data, 'Temperature', times[t_index_array], tube_id)
  dataManip.printSimpleStatsTubeData(F_data, 'Fuel', times[t_index_array], tube_id)
  dataManip.printSimpleStatsTubeData(Ma, 'MachNumber', times[t_index_array], tube_id)


  if 'PassiveScalars' in datas_dict:
    propTitles = ['temperature [K]', 'pressure [bar]', 'local Machnumber [-]', 'fuel Massfraction [-]', 'Passive Scalars [-]']
    titles = ['T [K]', 'P [bar]', 'Ma [-]', r'$\text{X}_{\text{fuel}}\,[-]$', 'Passive Scalars [-]']
    datas = [T_data * T_ref, p_data, Ma, F_data, X]
    f, ax = plt.subplots(nrows=1, ncols=5, figsize=(50. / 2, 10 / 2.), sharey=False)# figsize=(15, 10)) #set_size('thesis'))
  else:
    propTitles = ['temperature [K]', 'pressure [bar]', 'local Machnumber [-]', 'fuel Massfraction [-]']
    titles = [r'$T$'+' [K]', r'$p$'+' [bar]', r'$Ma$'+' [-]', r'$X_{fuel}$'+' [-]']
    datas = [T_data * T_ref, p_data, Ma, F_data]
    if CROPPED:
      f, ax = plt.subplots(nrows=2, ncols=2, figsize=(18. / 2., 5.), sharey=False)
    else:
      f, ax = plt.subplots(nrows=2, ncols=2, figsize=(18. / 2., 18. / 2.), sharey=False)
    
    ax = ax.flatten()
  def props(title):
    props = {
      'origin': 'lower',
      'interpolation': 'none',
      'extent': (x0, xEnd, t0, tEnd),
      'aspect': 'auto'
    }
    if title == 'Passive Scalars [-]':
      props = {
        'origin': 'lower',
        'extent': (x0, xEnd, t0, tEnd),
        'levels': np.linspace(np.min(X), np.max(X), 20),
        'vmax': None,
        'vmin': None,
        'cmap': 'twilight'
      }
    if title == 'fuel Massfraction [-]':
      props['vmax'] = 0.
      props['vmin'] = 1.
      props['cmap'] = 'gray_r'
    # if  title == 'temperature':
      # props['vmax'] = 3.0
    if title == 'pressure [bar]':
      props = {
      'origin': 'lower',
      'interpolation': 'none',
      'aspect': 'auto',
      'extent': (x0, xEnd, t0, tEnd),
      'vmin': 0.0,
      'vmax': 30.0,
      'cmap': 'jet',
      }
    return props
  import itertools
  ims = [a.imshow(data, **props(title)) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 4,  enumerate(zip(ax, datas, propTitles)))]
  # ims = [a.plot(data[-1,:]) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 4,  enumerate(zip(ax, datas, propTitles)))]
  if 'PassiveScalars' in datas_dict:
    ims.append(ax[4].contourf(datas[4], **props(propTitles[4])))
  for a, title in zip(ax, titles):
    a.set(ylabel=r'$t$'+' [-]')
    a.set(xlabel=r'$x$'+' [-]', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im, tit in zip(ax, ims, propTitles):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    a.label_outer()
    if 'pressure' in tit:
      plt.colorbar(im, ax=a, extend='max')
    else:
      plt.colorbar(im, ax=a)
  
  # if CROPPED:
  #   f.subplots_adjust(wspace=0.35, hspace=0.4)
  # else:
  #   f.suptitle("Tube id = {}".format(tube_id), y=0.93, fontsize=14)

  f.savefig(output_path+'/Tube{}-poster.png'.format(tube_id), bbox_inches='tight', dpi=175)
  f.savefig(output_path+'/Tube{}-poster.pdf'.format(tube_id), bbox_inches='tight', dpi=600)
  f.clear()
  plt.close(f)
  
  