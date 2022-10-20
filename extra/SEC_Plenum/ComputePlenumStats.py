import sys
import os
import numpy as np

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
import amrex.other as other
import amrex.h5_data_processing as dataManip

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# check cli
if len(sys.argv)<2:
  errMsg = ('Not enough input arguments!\n'
               +'\t1. argument must be dataPath!\n'
               +'\toptional argument is name of the inputfile\n'
               +'\toptional argument parallel working with module multiprocessing --parallel (default works on 8 CPUs)'
               +'\te.g. {} path --config=inputfile.py --parallel=4'.format(sys.argv[0])
            )
  raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

# name of the inputfile is optional
optional = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--config=' in el ]
if not optional:
    optional = ['inputfile.py'] # default value 
inputfileName = optional[0]

PARALLEL=False
if any('--parallel' in arg for arg in sys.argv):
   from multiprocessing import Pool, cpu_count
   PARALLEL=True
   Num_CPUs = 8
   try:
      cpuString = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--parallel' in el ]
      if cpuString:
         Num_CPUs = cpuString[0]
   except: 
      pass
   if Num_CPUs>cpu_count():
      Num_CPUs = cpu_count()
   print(f'starting computations on {Num_CPUs} from {cpu_count()} available cores')

other.import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import R #Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref

SYMMETRYCHECK = True
#-----------------------------------------------------------------------------

plenum = "{}/Plenum.h5".format(dataPath)
output_path = '{}/Visualization/'.format(dataPath)
os.makedirs(output_path, exist_ok=True)

print("Processing data in {}".format(output_path))

# Get times and nSteps from plenum
times = io.h5_load_timepoints(plenum)
nSteps = times.shape[0]

scalar_limits = io.getPassiveScalarLimits(plenum, 0, nSteps-1)

dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Pressure', 1.0, 8, output_path, True, PARALLEL)
dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Density', 1.0, 8, output_path, True, PARALLEL)
dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Temperature', 1.0, 8, output_path, True, PARALLEL)
dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'PassiveScalar', 1.0, 8, output_path, True, PARALLEL)

def getPlenumStats(i, time, PARALLEL=PARALLEL):
   if PARALLEL:
      print(f'sampling data from timepoint {round(time, 4)}')
   
   plenum_variables = ["Pressure", "Density", "PassiveScalars", 'vfrac']
   
   plenum_data, current_time, _, plenum_dict = io.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   
   volume_fraction = plenum_data[plenum_dict['vfrac']]
   pressure = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Pressure']], volume_fraction)
   rho = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Density']], volume_fraction)
   
   passiveScalarMF = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['PassiveScalars']], volume_fraction)
   
   temperature = pressure / (rho*R)

   # # print out the first occurence of min/max value 
   dataManip.printSimpleStatsPlenumSingleTimepoint(pressure, 'Pressure', current_time, output_path=output_path, PARALLEL=PARALLEL, SYMMETRYCHECK=SYMMETRYCHECK)
   dataManip.printSimpleStatsPlenumSingleTimepoint(rho, 'Density', current_time, output_path=output_path, PARALLEL=PARALLEL, SYMMETRYCHECK=SYMMETRYCHECK)
   dataManip.printSimpleStatsPlenumSingleTimepoint(temperature, 'Temperature', current_time, output_path=output_path, PARALLEL=PARALLEL, SYMMETRYCHECK=SYMMETRYCHECK)
   dataManip.printSimpleStatsPlenumSingleTimepoint(passiveScalarMF, 'PassiveScalar', current_time, output_path=output_path, PARALLEL=PARALLEL, SYMMETRYCHECK=SYMMETRYCHECK)

if PARALLEL:
   values = ((i, time, PARALLEL) for i, time in enumerate(times))

   with Pool(Num_CPUs) as pool:
      pool.starmap(getPlenumStats, values)
else:
   for i, timepoint in other.progressBar(times, enumeration=True):
      getPlenumStats(i, times)

# sort unorderd data and save them to disk
if PARALLEL:
   dataManip.sortSimpleStatsPlenum(output_path)
