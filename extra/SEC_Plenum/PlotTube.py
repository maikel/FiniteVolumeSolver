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
import h5py

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

tube_id = 1
#modes = ['cellwise', 'average_mirror_state', 'average_ghost_state']
#modes = ['cellwise']
modes = ['average_massflow']
#mode_titles = ['Cellwise', 'Average Mirror State', 'Average Ghost State']
mode_titles = ['Average Massflow']

for mode, suptitle in zip(modes, mode_titles):
  output_base_path = FVS_path+'/build_2D-Release/{}/Tube{}.h5'.format(mode, tube_id)
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
  T_data = p_data / rho_data
  F_data = rhoF_data / rho_data

  x0 = extent_1d[0]
  xEnd = extent_1d[1]
  t0 = 0.0
  tEnd = times[-1]

  # print(np.max(F_data))
  Ma = rhou_data / rho_data / c_data

  titles = ['Temperature', 'Pressure', 'Local Machnumber', 'Fuel Massfraction']
  datas = [T_data, p_data, Ma, F_data]
  f, ax = plt.subplots(nrows=1, ncols=4, figsize=(40. / 2, 10 / 2.), sharey=True)# figsize=(15, 10)) #set_size('thesis'))
  
  def props(title):
    props = {
      'origin': 'lower',
      'interpolation': 'none',
      'extent': (x0, xEnd, t0, tEnd),
      'aspect': 'auto'
    }
    if title == 'Fuel Massfraction':
      props['vmax'] = None
      props['vmin'] = None
    # if  title == 'Temperature':
      # props['vmax'] = 3.0
    if title == 'Pressure':
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

  ims = [a.imshow(data, **props(title)) for a, data, title in zip(ax, datas, titles)]
  ax[0].set(ylabel='time')
  for a, title in zip(ax, titles):
    a.set(xlabel='x', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im in zip(ax, ims):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.colorbar(im, ax=a)
  # f.suptitle(suptitle)
  f.savefig('Tube{}_{}.png'.format(tube_id, mode), bbox_inches='tight')
  f.clear()
  plt.close(f)
  