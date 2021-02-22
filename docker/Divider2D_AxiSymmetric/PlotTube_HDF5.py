import sys, os
# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
sys.path.append(pathname)
import amrex.plotfiles as da
#from amrex.plotfiles import h5_load_timeseries

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1])
   inputFilePath = dataPath
else:
   dataPath = 'ConvergentNozzle'

output_path = '{}/Visualization'.format(dataPath)

# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False

# inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

plt.style.use('seaborn')
# tex_fonts = {
#     # Use LaTeX to write all text
#     "text.usetex": True,
#     "font.family": "serif",
#     # Use 10pt font in plots, to match 10pt font in document
#     "axes.labelsize": 9,
#     "axes.titlesize": 9,
#     "axes.labelsize": 9,
#     "font.size": 9,
#     # Make the legend/label fonts a little smaller
#     "legend.fontsize": 9,
#     "xtick.labelsize": 7,
#     "ytick.labelsize": 7
# }

# plt.rcParams.update(tex_fonts)
plt.rcParams.update({'axes.grid' : False})

os.makedirs(output_path, exist_ok=True)

n_tubes = 1
for tube_id in range(n_tubes):
  print("Plotting Tube with id = {}".format(tube_id))
  filename_basic = dataPath+'/Tube.h5'
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

   print(fileNameList)

   # Read in data
   # Attention the last file is the latest!!
   # for example we have Filename.h5 and Filename.h5.1 the last one contains the first data!
   print("Read in data from {}".format(fileNameList[-1]))
   datas, times, datas_dict = da.h5_load_timeseries(fileNameList[-1])

   for filename in reversed(fileNameList[:-1]):
      print("Read in data from {}".format(filename))
      data, time, _ = da.h5_load_timeseries(fileNameList[0])
      datas = np.concatenate((datas, data))
      times = np.concatenate((times, time))
  else:
    print("Read in data from {}".format(filename_basic))
    datas, times, datas_dict = da.h5_load_timeseries(filename_basic)


  datas = np.squeeze(datas) # remove last axis
  print(datas.shape)
  # print(datas_dict)

  # slice_start = datas.shape[0]//3
  slice_start = 0
  
  T_ref = 1.0
  p_ref = 1.0

  rho = datas[slice_start:, datas_dict['Density'], :]
  rhou = datas[slice_start:, datas_dict['Momentum'], :]
  rhoE = datas[slice_start:, datas_dict['Energy'], :]
  rhoH2 = datas[slice_start:, datas_dict['H2'], :]
  p = datas[slice_start:, datas_dict['Pressure'], :]
  c = datas[slice_start:, datas_dict['SpeedOfSound'], :]
  rhoO2 = datas[slice_start:, datas_dict['O2'], :]
  T = datas[slice_start:, datas_dict['Temperature'], :]

  x0 = extent_1d[0]
  xEnd = extent_1d[1]
  t0 = times[slice_start]
  tEnd = times[-1]

  # print(np.max(F))
  u = rhou / rho
  O2 = rhoO2 / rho
  H2 = rhoH2 / rho

  titles = ['Temperature [K]', 'Pressure [bar]', 'Velocity [m/s]', 'H2 Massfraction [-]', 'O2 Massfraction [-]']
  datas = [T, p, u, H2, O2]
  f, ax = plt.subplots(nrows=1, ncols=5, figsize=(50. / 2, 10 / 2.), sharey=True)# figsize=(15, 10)) #set_size('thesis'))

  def props(title):
    props = {
      'origin': 'lower',
      'interpolation': 'none',
      'extent': (x0, xEnd, t0, tEnd),
      'aspect': 'auto',
      'cmap': 'turbo'
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
    if title == 'H2 Massfraction [-]' or title == 'O2 Massfraction [-]':
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
      # 'vmin': 0.7,
      # 'vmax': 1.4,
      'cmap': 'jet'
      }
    return props
  import itertools
  ims = [a.imshow(data, **props(title)) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 5,  enumerate(zip(ax, datas, titles)))]
  # ims = [a.plot(data[-1,:]) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 4,  enumerate(zip(ax, datas, titles)))]
  ax[0].set(ylabel='time')
  # ims.append(ax[4].contourf(datas[4], **props(titles[4])))
  for a, title in zip(ax, titles):
    a.set(xlabel='x', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im in zip(ax, ims):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # a.set(ylim=(150,None))
    plt.colorbar(im, ax=a)
  f.suptitle("Tube id = {}".format(tube_id))
  f.savefig(output_path+'/Tube{}.png'.format(tube_id), bbox_inches='tight')
  f.clear()
  plt.close(f)

  # f = plt.figure()
  # plt.plot(times, p_data[:,0], label='0')
  # plt.plot(times, p_data[:,1], label='1')
  # plt.plot(times, (p_data[:,0]+p_data[:,1])/2, label='mean')
  # # for i in range(10):
  # #   print("Cell {} pressure_max = {}".format(i, np.max(p_data[:,i])))

  
  # plt.ylim(5.85, 6.05)
  # plt.xlim(130,140)
  # plt.legend()
  # plt.grid(True)
  # f.savefig(output_path+'/Tube{}_pressureFirstCell.png'.format(tube_id), bbox_inches='tight')
  # f.clear()
  # plt.close(f)

  # f = plt.figure()
  # plt.plot(times, T_data[:,0], label='0')
  # plt.plot(times, T_data[:,1], label='1')
  # # plt.plot(times, (T_data[:,0]+T_data[:,1])/2, label='mean')
  # # plt.ylim(0, 1)
  # plt.xlim(130,140)
  # plt.legend()
  # plt.grid(True)
  # f.savefig(output_path+'/Tube{}_TemperatureFirstCell.png'.format(tube_id), bbox_inches='tight')
  # f.clear()
  
  