import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.plotfiles as da

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import h5py
import itertools

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1])
   inputFilePath = dataPath
else:
   dataPath = FVS_path+"/build_2D-Release/average_massflow"
   inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"

try:
  inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
  inputfileName = 'SEC_Plenum_Arrhenius.py'

da.import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref
from inputfile import D as diameter_tube

try:
  from inputfile import T_ref, n_tubes  
except:
  from inputfile import T_ref
  n_tubes=1

plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization/Plenum/'.format(outPath)

# Get times and nSteps from plenum
times = da.h5_load_timepoints(plenum)
nSteps = times.shape[0]

# optional slicing in time-dimension
tplotmin = 200.0
splittedPath = dataPath.rsplit('/',1)[-2] #-2 because / at the end
if 'vol40.0' in splittedPath:
   tplotmin = 200.0
elif 'vol20.0' in splittedPath:
   tplotmin = 145.0
tplotmax = 400.0
t_bool_array = (times>=tplotmin) & (times<=tplotmax)
t_index_array = np.nonzero(t_bool_array)[0] # returns tuple: (array([10 11 ...]), )

# check if index array is empty
if not np.any(t_index_array):
   raise IndexError('time index array is empty!')


# get data from output dictionary
output_dict = Output['outputs']
output_dict = [ out_dict for out_dict in output_dict if 'which_block' in out_dict] # strip 'CounterOutput'
for out_dict in output_dict:
   # get only the file names without the prefix dir names
   out_path = out_dict['path'].rsplit('/',1)[-1]
   if 'Plenum' in out_path:
      if 'frequencies' in out_dict:
         plenum_out = int(out_dict['frequencies'][0])
      elif 'intervals' in out_dict:
         plenum_out = float(out_dict['intervals'][0])
      else:
         raise Exception("could not find output frequency or interval from plenum")
   elif 'Tube' in out_path: # assume all tubes have same frequency or intervals
      if 'frequencies' in out_dict:
         tube_out = int(out_dict['frequencies'][0])
      elif 'intervals' in out_dict:
         tube_out = float(out_dict['intervals'][0])
      else:
         raise Exception("could not find output frequency or interval from tube")

# get the interval ratio from plenum and tubes
# normally we write the tube data 4 times more often than the plenum data
tube_output_factor = int(plenum_out / tube_out)
# print(tube_output_factor, plenum_out, tube_out)


def getPassiveScalarLimits(first=0, last=nSteps-1):
   (rho, rhoX, vols), _, _, _ = da.h5_load_spec_timepoint_variable(plenum, first, ["Density", "PassiveScalars", "vfrac"])
   rho = np.ma.masked_array(rho, vols < 1e-14)
   min = np.min(rhoX / rho)

   (rho, rhoX, vols), _, _, _ = da.h5_load_spec_timepoint_variable(plenum, last, ["Density", "PassiveScalars", "vfrac"])
   rho = np.ma.masked_array(rho, vols < 1e-14)
   max = np.max(rhoX / rho)

   return np.rint((min, max))

scalar_limits = getPassiveScalarLimits(first=t_index_array[0], last=t_index_array[-1])

def PrintProgress(i):
  progress = int(100.0 * float(i) / (nSteps - 1))
  print('[{:3d}%] Reading slice [{}/{}]'.format(progress, i + 1, nSteps))

os.makedirs(output_path, exist_ok=True)


da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Pressure', 1.0, output_path=output_path, firstCall=True)
da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Density', 1.0, output_path=output_path, firstCall=True)
da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Temperature', 1.0, output_path=output_path, firstCall=True)
da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'PassiveScalar', 1.0, output_path=output_path, firstCall=True)

for i in t_index_array:
   PrintProgress(i)
   
   plenum_variables = ["Pressure", "Density", "PassiveScalars", 'vfrac']
   
   plenum_data, current_time, plenum_extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   volume_fraction = plenum_data[plenum_dict['vfrac']]

   volume_fraction_offset = 1.e-14
   pressure = np.ma.masked_array(plenum_data[plenum_dict['Pressure']], volume_fraction < volume_fraction_offset)
   rho = np.ma.masked_array(plenum_data[plenum_dict['Density']], volume_fraction < volume_fraction_offset)
   passiveScalarMF = np.ma.masked_array(plenum_data[plenum_dict['PassiveScalars']], volume_fraction < volume_fraction_offset) / rho
   temperature = pressure / rho

   # # print out the first occurence of min/max value 
   da.printSimpleStatsPlenumSingleTimepoint(pressure, 'Pressure', current_time, output_path=output_path)
   da.printSimpleStatsPlenumSingleTimepoint(rho, 'Density', current_time, output_path=output_path)
   da.printSimpleStatsPlenumSingleTimepoint(temperature, 'Temperature', current_time, output_path=output_path)
   da.printSimpleStatsPlenumSingleTimepoint(passiveScalarMF, 'PassiveScalar', current_time, output_path=output_path)