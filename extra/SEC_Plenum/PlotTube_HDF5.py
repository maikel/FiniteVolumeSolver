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


inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"
sys.path.append(inputFilePath)
from SEC_Plenum_Arrhenius import T_ref, Tubes

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

n_tubes = len(Tubes)

dataPath = FVS_path+"/build_2D-ReleaseDebInfo/average_massflow"
output_path = '{}/Visualization'.format(dataPath)

for tube_id in range(n_tubes):
  print("Plotting Tube with id = {}".format(tube_id))
  output_base_path = dataPath+'/Tube{}.h5'.format(tube_id)
  datas, times, datas_dict = da.h5_load_timeseries(output_base_path)
  extent_1d = da.h5_load_get_extent_1D(output_base_path)
  
  datas = np.squeeze(datas) # remove last axis
  # print(datas.shape)
  # print(datas_dict)

  rho_data = datas[:, datas_dict['Density'], :]
  rhou_data = datas[:, datas_dict['Momentum'], :]
  rhoE_data = datas[:, datas_dict['Energy'], :]
  rhoF_data = datas[:, datas_dict['Species'], :]
  p_data = datas[:, datas_dict['Pressure'], :]
  c_data = datas[:, datas_dict['SpeedOfSound'], :]
  rhoX_data = datas[:, datas_dict['PassiveScalars'], :]
  T_data = p_data / rho_data
  F_data = rhoF_data / rho_data

  x0 = extent_1d[0]
  xEnd = extent_1d[1]
  t0 = 0.0
  tEnd = times[-1]

  # print(np.max(F_data))
  Ma = rhou_data / rho_data / c_data
  X = rhoX_data / rho_data

  titles = ['Temperature [K]', 'Pressure [bar]', 'Local Machnumber [-]', 'Fuel Massfraction [-]', 'Passive Scalars [-]']
  datas = [T_data * T_ref, p_data, Ma, F_data, X]
  f, ax = plt.subplots(nrows=1, ncols=5, figsize=(50. / 2, 10 / 2.), sharey=True)# figsize=(15, 10)) #set_size('thesis'))

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
      # 'vmin': 0.7,
      # 'vmax': 1.4,
      'cmap': 'jet'
      }
    return props
  import itertools
  ims = [a.imshow(data, **props(title)) for (__, (a, data, title)) in itertools.takewhile(lambda x: x[0] < 4,  enumerate(zip(ax, datas, titles)))]
  ax[0].set(ylabel='time')
  ims.append(ax[4].contourf(datas[4], **props(titles[4])))
  for a, title in zip(ax, titles):
    a.set(xlabel='x', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im in zip(ax, ims):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.colorbar(im, ax=a)
  f.suptitle("Tube id = {}".format(tube_id))
  f.savefig(output_path+'/Tube{}.png'.format(tube_id), bbox_inches='tight')
  f.clear()
  plt.close(f)
  