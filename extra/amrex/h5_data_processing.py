import numpy as np
import os
from datetime import datetime

def maskPlenumCutCells(data, vfrac, vfracOffset=1.0e-14):
  """
  Mask all Cut Cells from the (2D) Plenum data.

  Parameters
  ----------------------------------------
    data:         numpy array
                  data to mask
    vfrac:        numpy array (same shape as data)
                  volumefraction
    vfracOffset:  float
                  offset for mask

  Returns the masked numpy array
  """
  return np.ma.masked_array(data, vfrac < vfracOffset)

def printSimpleStatsTubeData(data, variable, times, tube_id=0, ndig=4, output_path=""):
  """
  Print out simple Stats from given Arrays. Only tested with HDF5 Tube Data!
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

  format_row = '{:>18}'*len(stats_data[0])
  print("[Tube{}] Stats for {}:".format(tube_id, variable))
  for row in stats_data:
    print(format_row.format(*row))
    if output_path:
      with open(fname, 'a') as f:
        f.write(format_row.format(*row)+"\n")
  print()

def printSimpleStatsPlenumSingleTimepoint(data, variable, time, ndig=8, output_path="", FIRSTCALL=False, PARALLEL=False):
  """
  Print out simple Stats from given Arrays.
  Shape must be (NCellsX, NCellsY, (NCellsZ))

  Parameters
  ----------------------------------------
    data:       numpy array
                the data array with shape like (NCellsX, NCellsY, (NCellsZ))
    variable:   string
                name of the variable
    times:      numpy array
                the time points from the data
    ndig:       integer
                maximal number of decimal digits for output
    output_path: string, optional
                 optional write out all stats in file
    FIRSTCALL:  bool, optional
                rename old file if any, write header in new file
  """
  if output_path:
    if 'Plenum' in output_path:
      fname = os.path.join( output_path.split('Plenum')[0], "Plenum_{}_stats.dat".format(variable) )
    else:
      fname = os.path.join( output_path, "Plenum_{}_stats.dat".format(variable) )
    if PARALLEL:
      fname = fname[:-4]+'_unorderd'+fname[-4:] # insert unorderd hint in fname!!
    if FIRSTCALL:
      if os.path.isfile(fname): # if file exists, rename it
        date = datetime.now().strftime("%Y_%m_%d-%I:%M:%S_%p")
        newName = fname.split('.dat')[0]+"_{}.dat".format(date)
        print("renamed old File '{}' to '{}'".format(fname, newName))
        os.rename(fname, newName)
      header = ['time', 'min', 'mean', 'median', 'std', 'max']
      format_row = '#{:>17}'+'{:>18}'*(len(header)-1)
      with open(fname, 'w') as f:
        f.write(format_row.format(*header)+"\n")
      return None
  
  if np.ma.is_masked(data):
    # version for masked numpy arrays
    indices_min = np.unravel_index(np.ma.argmin(data, axis=None), data.shape)
    indices_max = np.unravel_index(np.ma.argmax(data, axis=None), data.shape)
    stats_data = [time, data[indices_min], np.ma.mean(data), np.ma.median(data), np.ma.std(data), data[indices_max]]
  else: 
    indices_min = np.unravel_index(np.argmin(data, axis=None), data.shape)
    indices_max = np.unravel_index(np.argmax(data, axis=None), data.shape)
    stats_data = [time, data[indices_min], np.mean(data), np.median(data), np.std(data), data[indices_max]]
  
  stats_data = [el if isinstance(el, str) else round(el, ndig) for el in stats_data]
  
  format_row = '{:>18}'*len(stats_data)
  # print(format_row.format(*stats_data))
  if output_path:
    with open(fname, 'a') as f:
      f.write(format_row.format(*stats_data)+"\n")
  # print()