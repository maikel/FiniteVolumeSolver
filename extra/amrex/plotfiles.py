import yt
import numpy as np
import h5py

def yt_load(path_to_plotfile, vars, mask_boundary_cells=True, buffer_size=None, vfrac_cutoff=1e-15):
  """
  Load the with 'vars' specified 2D data from the plt-files.

  Parameters
  ----------------------------------------
    path_to_plotfile:     string
                          path to the plt-files
    vars:                 list of strings
                          list of variables to extract data
    mask_boundary_cells:  boolean, optional
                          mask the cutcell data with the volume fraction
    buffer_size:          optional
                          buffers loading data
    vfrac_cutoff:         float, optional
                          the cutoff parameter for masking the cutcell data

  Returns
  ----------------------------------------
    datas:        numpy array
                  the loaded datas from the plt-file
    current_time: float
                  the time point from the data
    bounds:       list of floats
                  the physical extension from the dataset
  """
  ds = yt.load(path_to_plotfile)
  plot = yt.plot_2d(ds, vars, origin='native')
  if buffer_size:
    plot.set_buff_size(buffer_size); 
  else:
    plot.set_buff_size(ds.domain_dimensions)
  current_time = float(ds.current_time)
  bounds = np.array(plot.bounds)
  datas = [np.squeeze(np.array(plot.frb[var])) for var in vars]
  if mask_boundary_cells:
    vfrac = np.array(plot.frb["vfrac"])
    vfrac = np.squeeze(vfrac)
    datas = [np.ma.masked_where( vfrac<=vfrac_cutoff, data) for data in datas]
  return datas, current_time, bounds

def h5_load_spec_timepoint_variable(filename, num, variables):
  """
  Load the with 'variables' specified 2D data from the HDF5 file for a given timepoint.

  Parameters
  ----------------------------------------
    filename:   string
                name of the HDF5-file
    num:        integer
                number of timepoint
    variables:  list of strings
                list of variables to extract data

  Returns
  ----------------------------------------
    datas:        tuple of numpy arrays
                  the loaded datas from the HDF5-file
    time:         float
                  the time point from the data
    extent:       list of floats
                  the physical extension from the dataset,
                  something like [xlo, xhi, ylo, yhi]
    dictionary:   dict
                  dictionary that matches the variable names to their index
  """
  file = h5py.File(filename, mode='r')
  strings = list(file['fields'].asstr())
  indices = [strings.index(var) for var in variables]
  dictionary = dict( zip(variables, range(len(variables))))
  shape = file['data'].shape
  nx = np.array([shape[2], shape[3]])
  datas = [np.reshape(np.array(file['data'][num, var, :, :]), [nx[1], nx[0]])  for var in indices]
  time = float(file['times'][num])
  dx = np.array(file['data'].attrs['cell_size'])
  xlower = np.array(file['data'].attrs['xlower'])
  xupper = xlower + dx * nx
  extent = [xlower[0], xupper[0], xlower[1], xupper[1]]
  file.close()
  return tuple(datas), time, extent, dictionary

def h5_load_timeseries(filename):
  """
  Load the all 1D data from the HDF5 file for all time points.

  Parameters
  ----------------------------------------
    filename:   string
                name of the HDF5-file
  Returns
  ----------------------------------------
    datas:        tuple of numpy arrays
                  the loaded datas from the HDF5-file
    time:         numpy array
                  the time points from the data
    dictionary:   dict
                  dictionary that matches the variable names to their index
  """
  file = h5py.File(filename, mode='r')
  strings = list(file['fields'].asstr())
  indices = [strings.index(var) for var in strings]
  dictionary = dict( zip(strings, indices) )
  shape = file['data'].shape
  datas = np.zeros((shape))
  for i ,var in enumerate(indices):
    datas[:,i] = file['data'][:, var]
  time = np.array(file['times'])
  file.close()
  return datas, time, dictionary

def h5_load_get_extent_1D(filename):
  """
  Load only the physical extension from the HDF5 file, assume it is an 1D dataset.

  Parameters
  ----------------------------------------
    filename:   string
                name of the HDF5-file

  Returns
  ----------------------------------------
    extent_1d:    list of floats
                  the physical extension from the 1D dataset,
                  something like [xlo, xhi]
  """
  file = h5py.File(filename, mode='r')
  shape = file['data'].shape
  nx = np.array([shape[2]])
  dx = np.array(file['data'].attrs['cell_size'])
  xlower = np.array(file['data'].attrs['xlower'])
  xupper = xlower + dx * nx
  extent_1d = [xlower[0], xupper[0]]
  file.close()
  return extent_1d

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