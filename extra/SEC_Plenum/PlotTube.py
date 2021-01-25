import numpy as np
import matplotlib 
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import h5py
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

def Load(path, chunk=[]):
    # @brief Reads grid data and time data from a specified hdf5 file.
    # @param path The path of the HDF5 file
    # @param dataset A string name of the main dataset containing the states data
    # @param times A string name of the time dataset containing all time points
    # @param chunk An optional [start, end] list that specifies how many time steps to load
    file = h5py.File(path, "r")
    if chunk:
        data_array = np.array(file['data'][chunk[0]:chunk[1],0,:,:])
        time_array = np.array(file['times'][chunk[0]:chunk[1]])
    else:
        data_array = np.array(file['data'][:, 0, :, :])
        time_array = np.array(file['times'])
    file.close()
    return data_array, time_array

  

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

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
  output_base_path = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Release/SEC_Plenum/{}/Tube{}.h5'.format(mode, tube_id)
  file = h5py.File(output_base_path, mode='r')
  rho_data = np.squeeze(np.array(file['data'][:, 0, :, :]))
  rhou_data = np.squeeze(np.array(file['data'][:, 1, :, :]))
  rhov_data = np.squeeze(np.array(file['data'][:, 2, :, :]))
  rhoE_data = np.squeeze(np.array(file['data'][:, 3, :, :]))
  rhoF_data = np.squeeze(np.array(file['data'][:, 4, :, :]))
  rhoFR_data = np.squeeze(np.array(file['data'][:, 5, :, :]))
  p_data = np.squeeze(np.array(file['data'][:, -2, :, :]))
  c_data = np.squeeze(np.array(file['data'][:, -1, :, :]))
  rhoP_data = rho_data - rhoF_data - rhoFR_data
  T_data = p_data / rho_data
  F_data = rhoF_data / rho_data
  times = np.squeeze(np.array(file['times']))
  file.close()

  x0 = -1.0
  xEnd = -0.09
  t0 = 0.0
  tEnd = times[-1]

  titles = ['Temperature', 'Pressure', 'Fuel Massfraction']
  datas = [T_data, p_data, F_data]
  f, ax = plt.subplots(nrows=1, ncols=3, figsize=set_size('thesis'))# figsize=(15, 10)) #set_size('thesis'))
  
  def props(title):
    props = {
      'origin': 'lower',
      'interpolation': 'none',
      'extent': (x0, xEnd, t0, tEnd),
      'aspect': 'auto'
    }
    if  title == 'Temperature':
      props['vmax'] = 3.0
    if title == 'Pressure':
      props = {
      'origin': 'lower',
      'interpolation': 'none',
      'aspect': 'auto',
      'extent': (x0, xEnd, t0, tEnd),
      'vmin': 0.7,
      'vmax': 1.4,
      'cmap': 'jet'
      }
    return props

  ims = [ax[0].imshow(datas[0], **props(titles[0])),
         ax[1].imshow(datas[1], **props(titles[1])),
         ax[2].imshow(datas[2], **props(titles[2]))]
  ax[0].set(ylabel='time')
  for a, title in zip(ax, titles):
    a.set(xlabel='x', title=title)

  from matplotlib.ticker import FormatStrFormatter
  for a, im in zip(ax, ims):
    a.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.colorbar(im, ax=a)
  f.suptitle(suptitle)
  f.savefig('Tube{}_{}.png'.format(tube_id, mode))
  f.clear()
  plt.close(f)
  