import sys, os
# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
sys.path.append(FVS_path+'/extra/')
import amrex.plotfiles as da
#from amrex.plotfiles import h5_load_timeseries

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1]) # path to data
   inputFilePath = dataPath # assumes inputfile is located in datapath
else:
   dataPath = FVS_path+"/build_2D-Release/average_massflow"
   inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"

output_path = '{}/Visualization'.format(dataPath)

# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False

try:
  inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
  inputfileName = 'SEC_Plenum_Arrhenius.py'

da.import_file_as_module(inputFilePath+inputfileName, 'inputfile')
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

os.makedirs(output_path, exist_ok=True)


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
  tplotmin = 0.0
  tplotmax = 400.0
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
  da.printSimpleStatsTubeData(p_data, 'Pressure', times[t_index_array], tube_id)
  da.printSimpleStatsTubeData(T_data, 'Temperature', times[t_index_array], tube_id)
  da.printSimpleStatsTubeData(F_data, 'Fuel', times[t_index_array], tube_id)
  da.printSimpleStatsTubeData(Ma, 'MachNumber', times[t_index_array], tube_id)


  if 'PassiveScalars' in datas_dict:
    titles = ['Temperature [K]', 'Pressure [bar]', 'Local Machnumber [-]', 'Fuel Massfraction [-]', 'Passive Scalars [-]']
    datas = [T_data * T_ref, p_data, Ma, F_data, X]
    f, ax = plt.subplots(nrows=1, ncols=5, figsize=(50. / 2, 10 / 2.), sharey=False)# figsize=(15, 10)) #set_size('thesis'))
  else:
    titles = ['Temperature [K]', 'Pressure [bar]', 'Local Machnumber [-]', 'Fuel Massfraction [-]']
    datas = [T_data * T_ref, p_data, Ma, F_data]
    f, ax = plt.subplots(nrows=1, ncols=4, figsize=(40. / 2, 10 / 2.), sharey=False)# figsize=(15, 10)) #set_size('thesis'))

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
    if title == 'Fuel Massfraction [-]':
      props['vmax'] = None
      props['vmin'] = None
      props['cmap'] = 'gray_r'
    # if  title == 'Temperature':
      # props['vmax'] = 3.0
    if title == 'Pressure [bar]':
      props = {
      'origin': 'lower',
      'interpolation': 'none',
      'aspect': 'auto',
      'extent': (x0, xEnd, t0, tEnd),
      'vmin': 0.0,
      'vmax': 30.0,
      'cmap': 'jet'
      }
    return props
  import itertools
  ims = [a.imshow(data, **props(title)) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 4,  enumerate(zip(ax, datas, titles)))]
  # ims = [a.plot(data[-1,:]) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 4,  enumerate(zip(ax, datas, titles)))]
  if 'PassiveScalars' in datas_dict:
    ims.append(ax[4].contourf(datas[4], **props(titles[4])))
  for a, title in zip(ax, titles):
    a.set(xlabel='x', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im in zip(ax, ims):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.colorbar(im, ax=a)
  f.suptitle("Tube id = {}".format(tube_id))
  f.savefig(output_path+'/Tube{}.png'.format(tube_id), bbox_inches='tight')
  f.savefig(output_path+'/Tube{}.pdf'.format(tube_id), bbox_inches='tight')
  f.clear()
  plt.close(f)

  # f = plt.figure()
  # plt.plot(times[:slice_end], p_data[:,0], label='0')
  # plt.plot(times[:slice_end], p_data[:,1], label='1')
  # # plt.plot(times, (p_data[:,0]+p_data[:,1])/2, label='mean')
  # # for i in range(10):
  # #   print("Cell {} pressure_max = {}".format(i, np.max(p_data[:,i])))

  
  # plt.ylim(5.85, 6.05)
  # plt.xlim(135.5,None)
  # plt.legend()
  # plt.grid(True)
  # f.savefig(output_path+'/Tube{}_pressureFirstCell.png'.format(tube_id), bbox_inches='tight')
  # f.clear()
  # plt.close(f)

  # f = plt.figure()
  # plt.plot(times[:slice_end], F_data[:,0], label='0')
  # plt.plot(times[:slice_end], F_data[:,1], label='1')
  # # plt.plot(times, (T_data[:,0]+T_data[:,1])/2, label='mean')
  # # plt.ylim(0, 1)
  # plt.xlim(135.5,None)
  # plt.legend()
  # plt.grid(True)
  # f.savefig(output_path+'/Tube{}_FuelFirstCell.png'.format(tube_id), bbox_inches='tight')
  # f.clear()
  
  