import sys
import os
import numpy as np

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as da
import amrex.other as other
import amrex.h5_data_processing as dataManip

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

other.import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref, R
from inputfile import D as diameter_tube

#-----------------------------------------------------------------------------

plenum = "{}/Plenum.h5".format(dataPath)
output_path = '{}/Visualization/'.format(dataPath)
os.makedirs(output_path, exist_ok=True)

print("Processing data in {}".format(output_path))

# Get times and nSteps from plenum
times = da.h5_load_timepoints(plenum)
nSteps = times.shape[0]

def getPassiveScalarLimits(plenumFile, first=0, last=nSteps-1):
   (rho, rhoX, vols), _, _, _ = da.h5_load_spec_timepoint_variable(plenumFile, first, ["Density", "PassiveScalars", "vfrac"])
   rho = np.ma.masked_array(rho, vols < 1e-14)
   min = np.min(rhoX / rho)

   (rho, rhoX, vols), _, _, _ = da.h5_load_spec_timepoint_variable(plenumFile, last, ["Density", "PassiveScalars", "vfrac"])
   rho = np.ma.masked_array(rho, vols < 1e-14)
   max = np.max(rhoX / rho)

   return np.rint((min, max))

scalar_limits = getPassiveScalarLimits(plenum)

dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Pressure', 1.0, output_path=output_path, firstCall=True)
dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Density', 1.0, output_path=output_path, firstCall=True)
dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Temperature', 1.0, output_path=output_path, firstCall=True)
dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'PassiveScalar', 1.0, output_path=output_path, firstCall=True)

for i, timepoint in other.progressBar(times, enumeration=True):

   plenum_variables = ["Pressure", "Density", "PassiveScalars", 'vfrac']
   
   plenum_data, current_time, plenum_extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   
   volume_fraction = plenum_data[plenum_dict['vfrac']]
   pressure = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Pressure']], volume_fraction)
   rho = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['Density']], volume_fraction)
   
   passiveScalarMF = dataManip.maskPlenumCutCells(plenum_data[plenum_dict['PassiveScalars']], volume_fraction)
   
   temperature = pressure / (rho*R)

   # # print out the first occurence of min/max value 
   dataManip.printSimpleStatsPlenumSingleTimepoint(pressure, 'Pressure', current_time, output_path=output_path)
   dataManip.printSimpleStatsPlenumSingleTimepoint(rho, 'Density', current_time, output_path=output_path)
   dataManip.printSimpleStatsPlenumSingleTimepoint(temperature, 'Temperature', current_time, output_path=output_path)
   dataManip.printSimpleStatsPlenumSingleTimepoint(passiveScalarMF, 'PassiveScalar', current_time, output_path=output_path)