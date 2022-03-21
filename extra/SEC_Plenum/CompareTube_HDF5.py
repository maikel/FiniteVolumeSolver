import sys, os
# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as da
import amrex.h5_data_processing as dataManip
import amrex.other as other

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import glob

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

other.import_file_as_module(os.path.join(inputFilePath, inputfileName), 'inputfile')
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
  # FVS_data, FVS_times, FVS_dict = da.h5_load_timeseries(filename_basic)
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
   FVS_data, FVS_times, FVS_dict = da.h5_load_timeseries(fileNameList[-1])

   for filename in reversed(fileNameList[:-1]):
      print("[Tube{}] Read in data from {}".format(tube_id, filename))
      data, time, _ = da.h5_load_timeseries(fileNameList[0])
      FVS_data = np.concatenate((FVS_data, data))
      FVS_times = np.concatenate((FVS_times, time))
  else:
    print("[Tube{}] Read in data from {}".format(tube_id, filename_basic))
    FVS_data, FVS_times, FVS_dict = da.h5_load_timeseries(filename_basic)


  FVS_data = np.squeeze(FVS_data) # remove last axis
  print("[Tube{}] data shape from tube is {} = (NTimePoints, NVariables, NCells)".format(tube_id, FVS_data.shape))
  # print(FVS_dict)


  ###### Get SEC_C data
  ## all possible directories with converted HDF5 files
  # hdf5_dir_full_list = ["p", "qc", "qr", "qv", "rho", "rhoe", "rhou", "rhov", "rhow", "rhoY", "S", "T", "Turb_Eps", "Turb_K", "u", "v", "w", "Y", "rhoX"]
  hdf5_dir_list = ["p", "rho", "rhou", "rhoX"]
  # hdf5_dir_list = ["Y"]
  hdf5_extension = '.h5'

  GT_dict = da.h5_get_Dict_Klein(hdf5_dir_list)

  GT_data=[]

  for dir in hdf5_dir_list:
    print("read in all HDF5 data from {}/{}".format(dataPath, dir))
    file_list = glob.glob("{}/{}/*{}".format(dataPath, dir, hdf5_extension))
    if not file_list:
      raise Exception("File list is empty! Check if '{}' is the right directory name".format(dir))
    file_list.sort()
    # print(file_list)

    temp_data = []
    GT_times = []
    
    for file in other.progressBar(file_list):
      # print("reading file {}".format(file))
      data, time, varname = da.h5_load_timeseries_Klein(file)

      temp_data.append(data)
      GT_times.append(time)
    GT_data.append(temp_data)

  GT_data = np.array(GT_data)
  GT_times = np.array(GT_times)
  GT_times = GT_times[1:-1]
  GT_data = np.swapaxes(GT_data,0,1)
  GT_data = GT_data[1:-1,:,2:-2]
  print("[Tube{} GT] data shape from tube is {} = (NTimePoints, NVariables, NCells)".format(tube_id, GT_data.shape))
  # print(GT_data.shape)
  gamma=1.4

  tplotmin = 100.0
  tplotmax = 400.0
  print("timeslicing ACTIVE tplotmin = {} and tplotmax={}".format(tplotmin, tplotmax))
  GT_t_index = (GT_times>=tplotmin) & (GT_times<=tplotmax)
  # print(GT_times[-5:])
  print("GT t0= {}; tend = {}".format(GT_times[GT_t_index][0], GT_times[GT_t_index][-1]))
  # 2:-2 --> skip ghost cell data!
  GT_rho = GT_data[GT_t_index, GT_dict['Density']]
  GT_rhou = GT_data[GT_t_index, GT_dict['Momentum_0']]
  GT_rhoF = GT_data[GT_t_index, GT_dict['Species']]
  GT_pressure = GT_data[GT_t_index, GT_dict['Pressure']]
  GT_c = np.sqrt( gamma * GT_pressure / GT_rho)
  GT_temperature = GT_pressure / GT_rho
  GT_fuel = GT_rhoF / GT_rho
  GT_mach = GT_rhou / GT_rho / GT_c
  ###### Get SEC_C data end



  # optional slicing in time-dimension
  FVS_t_index = (FVS_times>=tplotmin) & (FVS_times<=tplotmax)
  print("FVS t0= {}; tend = {}".format(FVS_times[FVS_t_index][0], FVS_times[FVS_t_index][-1]))
  # FVS_t_index = (FVS_times>=tplotmin)
  
  FVS_rho = FVS_data[FVS_t_index, FVS_dict['Density'], :]
  FVS_rhou = FVS_data[FVS_t_index, FVS_dict['Momentum'], :]
  FVS_rhoE = FVS_data[FVS_t_index, FVS_dict['Energy'], :]
  FVS_rhoF = FVS_data[FVS_t_index, FVS_dict['Species'], :]
  FVS_pressure = FVS_data[FVS_t_index, FVS_dict['Pressure'], :]
  FVS_c = FVS_data[FVS_t_index, FVS_dict['SpeedOfSound'], :]
  FVS_temperature = FVS_pressure / FVS_rho
  FVS_fuel = FVS_rhoF / FVS_rho
  FVS_mach = FVS_rhou / FVS_rho / FVS_c

  x0 = extent_1d[0]
  xEnd = extent_1d[1]
  t0 = FVS_times[FVS_t_index][0]
  tEnd = round( FVS_times[FVS_t_index][-1], 2) # round tEnd to display the last yticklabel
  print("[Tube{}] tEnd is {}".format(tube_id, tEnd))

  # # print out the first occurence of min/max value 
  # dataManip.printSimpleStatsTubeData(FVS_pressure, 'Pressure', FVS_times[t_index_array], tube_id)
  # dataManip.printSimpleStatsTubeData(FVS_temperature, 'Temperature', FVS_times[t_index_array], tube_id)
  # dataManip.printSimpleStatsTubeData(FVS_fuel, 'Fuel', FVS_times[t_index_array], tube_id)
  # dataManip.printSimpleStatsTubeData(FVS_mach, 'MachNumber', FVS_times[t_index_array], tube_id)
  
  def multi_interp(x, xp, fp):
    new_f = [ np.interp(x, xp, fp[:,i]) for i in range(fp.shape[1]) ]
    new_f = np.array(new_f)
    new_f = np.swapaxes(new_f, 0, 1)
    return new_f

  # GT_temperature = multi_interp(FVS_times[FVS_t_index], GT_times[GT_t_index], GT_temperature)
  # GT_pressure = multi_interp(FVS_times[FVS_t_index], GT_times[GT_t_index], GT_pressure)
  # GT_mach = multi_interp(FVS_times[FVS_t_index], GT_times[GT_t_index], GT_mach)
  # GT_fuel = multi_interp(FVS_times[FVS_t_index], GT_times[GT_t_index], GT_fuel)

  # FVS_temperature = multi_interp(GT_times[GT_t_index], FVS_times[FVS_t_index], FVS_temperature)
  # FVS_pressure = multi_interp(GT_times[GT_t_index], FVS_times[FVS_t_index], FVS_pressure)
  # FVS_mach = multi_interp(GT_times[GT_t_index], FVS_times[FVS_t_index], FVS_mach)
  # FVS_fuel = multi_interp(GT_times[GT_t_index], FVS_times[FVS_t_index], FVS_fuel)
  
  def flatten_a_python_list(list):
      # only for list in list [ [bla], [blab] ]
      return [item for sublist in list for item in sublist]

  def props(title):
    props = {
      'origin': 'lower',
      'interpolation': 'none',
      'extent': (x0, xEnd, t0, tEnd),
      'aspect': 'auto'
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
      'vmax': 10.0,
      'cmap': 'jet'
      }
    return props

  titles = ['Temperature [K]', 'Pressure [bar]', 'Local Machnumber [-]', 'Fuel Massfraction [-]']
  # Plot_data = [np.abs(FVS_temperature-GT_temperature), 
  #           np.abs(FVS_pressure-GT_pressure), 
  #           np.abs(FVS_mach-GT_mach), 
  #           np.abs(FVS_fuel-GT_fuel)]

  new_titles = []
  # new_titles.append(['FVS-GT: '+tit for tit in titles])
  new_titles.append(['FVS: '+tit for tit in titles])
  new_titles.append(['SEC\_C: '+tit for tit in titles])
  nrows=len(new_titles)
  
  titles = flatten_a_python_list(new_titles)
  
  new_data = []
  # new_data.append(Plot_data)
  new_data.append([FVS_temperature*T_ref, FVS_pressure, FVS_mach, FVS_fuel])
  new_data.append([GT_temperature*T_ref, GT_pressure, GT_mach, GT_fuel])
  Plot_data = flatten_a_python_list(new_data)
  # print(len(Plot_data))
  # print(Plot_data[5])

  
  f, ax = plt.subplots(nrows=nrows, ncols=4, figsize=(40. / 2, nrows*10*1.5 / 2.), sharey=False)# figsize=(15, 10)) #set_size('thesis'))

  import itertools
  ims = [a.imshow(data, **props(title.split(' ',1)[-1])) for (__, (a, data, title)) in enumerate(zip(ax.flatten(), Plot_data, titles))]
  
  for a, title in zip(ax.flatten(), titles):
    a.set(xlabel='x', ylabel='time', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im in zip(ax.flatten(), ims):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.colorbar(im, ax=a)
  # f.suptitle("Tube id = {}".format(tube_id))
  f.savefig(output_path+'/Tube_Compare_t0_{}_tend_{}.pdf'.format(int(t0), int(tEnd)), bbox_inches='tight')
  f.clear()
  plt.close(f)

  
  