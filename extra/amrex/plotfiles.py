import numpy as np
import h5py
import os

def yt_load(path_to_plotfile, vars, mask_boundary_cells=True, buffer_size=None, vfrac_cutoff=1e-15):
  import yt
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
  file = h5py.File(filename, mode='r', swmr=True)
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
  file = h5py.File(filename, mode='r', swmr=True)
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
  file = h5py.File(filename, mode='r', swmr=True)
  shape = file['data'].shape
  nx = np.array([shape[2]])
  dx = np.array(file['data'].attrs['cell_size'])
  xlower = np.array(file['data'].attrs['xlower'])
  xupper = xlower + dx * nx
  extent_1d = [xlower[0], xupper[0]]
  file.close()
  return extent_1d

def progressBar(iterable, prefix = 'Progress:', suffix = 'Complete', decimals = 1, length = 50, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    This fancy function is taken from https://stackoverflow.com/a/34325723.

    Parameters:
    ------------------------------
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\\r", "\\r\\n") (Str)

    Usage example
    ----------------------------
    \# create an iterable list of items \n
    items = list(range(0, 57))

    for item in progressBar(items, prefix = 'Progress:', suffix = 'Complete', length = 50):
    \n\# Do stuff...
    
    """
    total = len(iterable)
    # Progress Bar Printing Function
    def printProgressBar (iteration):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Initial Call
    printProgressBar(0)
    # Update Progress Bar
    for i, item in enumerate(iterable):
        yield item
        printProgressBar(i + 1)
    # Print New Line on Complete
    print("\n")

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

def get_controlState_Klein(filename):
  """
  Load the Control State data from the txt file for all time points.

  Parameters
  ----------------------------------------
    filename:   string
                name of the txt-file, usually GT_time_series.txt
  Returns
  ----------------------------------------
    GT_data:      numpy array
                  the loaded data from the txt-file
    times:        numpy array
                  the time points from the data
    dictionary:   dictionary
                  dict that matches the name to the data
  """
  ## all headers from GT
  # GT_header_list = ["time", "rpm", "PowV", "PowT", "PowOut", "ftot", "frate", "Eff", "MFVout", "MFVin", "MFTin", "MFTout", "rhoV", "pV", "TV", "rhoT", "pT", "TT", "pTmax"]
  # tcorresponding names for FVS Simulation
  FVS_header = ["time", "current_rpm", "compressor_power", "turbine_power", "power_out", "fuel_consumption", "fuel_consumption_rate", "efficiency", "compressor_mass_flow_out", "compressor_mass_flow_in", "turbine_mass_flow_in", "turbine_mass_flow_out", "compressor_density", "compressor_pressure", "compressor_temperature", "turbine_density", "turbine_pressure", "turbine_temperature", "pTmax"]
  
  # create dict without "time"
  dictionary = {key: value for value, key in enumerate(FVS_header[1:]) }

  GT_data = np.loadtxt(filename, skiprows=1) # first row contains header
  times = GT_data[:,0]
  GT_data = GT_data[:,1:]

  return GT_data, times, dictionary

def h5_get_Dict_Klein(hdf5_dir_list):
  """
  Create a dictionary corresponding to the given directory names.

  Parameters
  ----------------------------------------
    hdf5_dir_list:   list
                     list of directories to translate to dict
  Returns
  ----------------------------------------
    dictionary:      dictionary
                     the created dictionary
  """
  # all possible directory names from the simulation 
  hdf5_dir_full_list = ["p", "qc", "qr", "qv", "rho", "rhoe", "rhou", "rhov", "rhow", "rhoY", "S", "T", "Turb_Eps", "Turb_K", "u", "v", "w", "Y", "rhoX"]
  
  # the translated dictionary keys corresponding to the  directory names
  dict_list = ["Pressure", "qc", "qr", "qv", "Density", "Energy", "Momentum_0", "Momentum_1", "Momentum_2", "RhoTheta", "S", "Temperature", "Turb_Eps", "Turb_K", "Velocity_0", "Velocity_1", "Velocity_2", "Theta", "Species"]
  
  # find the indices from the directory list
  index_list = [hdf5_dir_full_list.index(element) for element in hdf5_dir_list if element in hdf5_dir_full_list]
  
  # and translate them to the dictionary
  dictionary = {dict_list[ind]: i for i, ind in enumerate(index_list) }

  return dictionary

def h5_load_timeseries_Klein(filename):
  """
  Load the all 1D data from the HDF5 file for all time points.

  Parameters
  ----------------------------------------
    filename:   string
                name of the HDF5-file
  Returns
  ----------------------------------------
    data:         numpy array
                  the loaded data from the HDF5-file
    time:         float
                  the time point from the data
    varname:      string
                  the name of the variable
  """
  file = h5py.File(filename, mode='r', swmr=True)
  data = np.array(file['Data-Set-2'])
  time = float(file['Data-Set-2'].attrs['valid_max'])
  varname = file['Data-Set-2'].attrs['long_name'].decode('UTF-8') # must be decoded because this is a byte string
  file.close()
  return data, time, varname

def printSimpleStatsTubeData(data, variable, times, tube_id=0, ndig=4, output_path=""):
  """
  Print out simple Stats from given Arrays. Only tested with Tube Data!
  Shape must be (NTimePoints, NCells)

  Parameters
  ----------------------------------------
    data:       numpy array
                the data array with shape like (NTimePoints, NCells)
    variable:   string
                name of the variable
    times:      numpy array
                the time points from the data
    tube_id:    integer, optional
                number of the tube, default=0
    ndig:       integer
                maximal number of decimal digits for output
    output_path: string, optional
                 optional write out all stats in file
  """
  if output_path:
    fname = os.path.join( output_path, "Tube{}_stats.dat".format(tube_id) )
    if os.path.isfile(fname):
      print("appending Tube Stats to File: {}".format(fname))
    with open(fname, 'a') as f:
      f.write("[Tube{}] Stats for {}:\n".format(tube_id, variable))

  indices_min = np.unravel_index(np.argmin(data, axis=None), data.shape)
  indices_max = np.unravel_index(np.argmax(data, axis=None), data.shape)

  stats_data = [['', 'min', 'mean', 'median', 'std', 'max'],
                ['value', data[indices_min], np.mean(data), np.median(data), np.std(data), data[indices_max]],
                ['time', times[indices_min[0]], '-', '-', '-', times[indices_max[0]]] ]
  
  stats_data[1] = [el if isinstance(el, str) else round(el, ndig) for el in stats_data[1]]
  stats_data[2] = [el if isinstance(el, str) else round(el, ndig) for el in stats_data[2]]

  format_row = '{:>12}'*len(stats_data[0])
  print("[Tube{}] Stats for {}:".format(tube_id, variable))
  for row in stats_data:
    print(format_row.format(*row))
    if output_path:
      with open(fname, 'a') as f:
        f.write(format_row.format(*row)+"\n")
  print()

def import_file_as_module(module_path, module_name):
  """
    Import a regular file as python module, following the python syntax: \n
    import module_path as module_name

    Parameters:
    ------------------------------
        module_path   - Required  : path to file with filename (str)
        module_name   - Required  : name of the module (str)

    Usage example
    ----------------------------
    \# import the inputfile as module inputfile \n
    import_file_as_module(inputFilePath+'SEC_Plenum_Arrhenius.py', 'inputfile') \n
    \# now we can import some objects from this new module \n
    from inputfile import T_ref, ControlOptions \n
    print(ControlOptions)
  """
  import sys, os
  if not os.path.isfile(module_path):
   raise FileNotFoundError('given file: {} does not exist!'.format(module_path))

  import importlib.machinery
  import importlib.util
  loader = importlib.machinery.SourceFileLoader(module_name, module_path)
  spec = importlib.util.spec_from_loader(loader.name, loader)
  mod = importlib.util.module_from_spec(spec)
  loader.exec_module(mod)
  sys.modules[module_name] = mod
  return mod