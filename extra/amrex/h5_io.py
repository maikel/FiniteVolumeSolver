import numpy as np
import h5py

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
    datas:        numpy array
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

def h5_load_timepoints(filename):
  """
  Load the all 1D data from the HDF5 file for all time points.

  Parameters
  ----------------------------------------
    filename:   string
                name of the HDF5-file
  Returns
  ----------------------------------------
    time:         numpy array
                  the time points from the data
  """
  file = h5py.File(filename, mode='r', swmr=True)
  time = np.array(file['times'])
  file.close()
  return time

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
