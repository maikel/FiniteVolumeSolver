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

# check cli
if len(sys.argv)<2:
   errMsg = ('Not enough input arguments!\n'
               +'\tfirst argument must be dataPath!')
   raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

try:
  inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
  inputfileName = 'SEC_Plenum_Arrhenius.py'

da.import_file_as_module(os.path.join(inputFilePath, inputfileName), 'inputfile')
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

titles = ['Pressure', 'Temperature',  'Passive Scalar', 'Velocity']
def props(title):
   nlevels = 96
   props = {
      # 'origin': 'lower',
      # 'interpolation': 'none',
      # 'aspect': 'equal',
      'extend': 'both',
   }
   if title == 'Passive Scalar':
      levels = np.linspace(*scalar_limits, nlevels)
      props.update(
         {
         'levels': levels,
         'vmin': levels[0],
         'vmax': levels[-1],
         'cmap': 'twilight'
         }
      )
   if  title == 'Temperature':
      minT = 5.*T_ref
      maxT = 12.5*T_ref
      scaleT = 25
      levels = np.linspace(minT, maxT, int(maxT-minT)//scaleT )
      props.update(
         {
            'levels': levels,
            'vmin': levels[0],
            'vmax': levels[-1],
            'cmap': 'afmhot_r'
         }
      )
      # props = {
      #    'origin': 'lower',
      #    'interpolation': 'none',
      #    'aspect': 'auto',
      #    # 'levels': np.linspace(10.0*T_ref, 13.0*T_ref, 30),
      #    'vmin': 6.5*T_ref,
      #    'vmax': 10.*T_ref,
      #    # 'cmap': 'jet'
      #    }
   if title == 'Pressure':
      levels = np.linspace(1.0, 20., nlevels)
      props.update(
         {
         'levels': levels,
         'vmin': levels[0],
         'vmax': levels[-1],
         'cmap': 'jet'
         }
      )
   if title == 'Velocity':
      props = {
         'angles': 'xy', 
         'scale_units': 'xy', 
         'scale': 10.0
      }
   return props


def PrintProgress(i):
  progress = int(100.0 * float(i) / (nSteps - 1))
  print('[{:3d}%] Reading slice [{}/{}]'.format(progress, i + 1, nSteps))

os.makedirs(output_path, exist_ok=True)

tube_paths = ['{}/Tube{}.h5'.format(dataPath, tube) for tube in range(n_tubes)]

def stackTubeDataTo2D(Tube_datalist):
   # all Tubedata is 1D but for contourf we need at least 2D data. so simply stack twice the 1d array
   for i, el in enumerate(Tube_datalist):
      el = np.squeeze(el)
      Tube_datalist[i] = np.stack((el,el))
   return Tube_datalist

# for i in itertools.dropwhile(lambda x: x < 53, range(nSteps)):
# for i in itertools.dropwhile(lambda i: i < 4997, range(nSteps)):

da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Pressure', 1.0, output_path=output_path, firstCall=True)
da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Density', 1.0, output_path=output_path, firstCall=True)
da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Temperature', 1.0, output_path=output_path, firstCall=True)
da.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'PassiveScalar', 1.0, output_path=output_path, firstCall=True)

for i in t_index_array:
   PrintProgress(i)
   
   Tube_p = []
   Tube_X = []
   tube_extents = []
   tube_variables = ["Pressure", "PassiveScalars", "Density"]

   for tube in tube_paths:
      tube_data, current_time, ext, tube_dict = da.h5_load_spec_timepoint_variable(tube, tube_output_factor*i, tube_variables)
      Tube_p.append(tube_data[tube_dict['Pressure']])
      Tube_X.append( tube_data[tube_dict['PassiveScalars']] / tube_data[tube_dict['Density']] )
      tube_extents.append(ext)
   
   Tube_p = stackTubeDataTo2D(Tube_p)
   Tube_X = stackTubeDataTo2D(Tube_X)
   
   plenum_variables = ["Pressure", "Density", "Momentum_0", "Momentum_1", "PassiveScalars", 'vfrac']
   # # alternative:
   # (pressure, rho, rhoU, rhoV, rhoX, vols), current_time, plenum_extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   plenum_data, current_time, plenum_extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   volume_fraction = plenum_data[plenum_dict['vfrac']]

   volume_fraction_offset = 1.e-14
   pressure = np.ma.masked_array(plenum_data[plenum_dict['Pressure']], volume_fraction < volume_fraction_offset)
   rho = np.ma.masked_array(plenum_data[plenum_dict['Density']], volume_fraction < volume_fraction_offset)
   rhoU = np.ma.masked_array(plenum_data[plenum_dict['Momentum_0']], volume_fraction < volume_fraction_offset)
   rhoV = np.ma.masked_array(plenum_data[plenum_dict['Momentum_1']], volume_fraction < volume_fraction_offset)
   passiveScalarMF = np.ma.masked_array(plenum_data[plenum_dict['PassiveScalars']], volume_fraction < volume_fraction_offset) / rho
   temperature = pressure / rho

   # # print out the first occurence of min/max value 
   da.printSimpleStatsPlenumSingleTimepoint(pressure, 'Pressure', current_time, output_path=output_path)
   da.printSimpleStatsPlenumSingleTimepoint(rho, 'Density', current_time, output_path=output_path)
   da.printSimpleStatsPlenumSingleTimepoint(temperature, 'Temperature', current_time, output_path=output_path)
   da.printSimpleStatsPlenumSingleTimepoint(passiveScalarMF, 'PassiveScalar', current_time, output_path=output_path)


   f, axs = plt.subplots(nrows=2, ncols=2, figsize=(20. / 2, 15. / 2), gridspec_kw={'width_ratios': [4,2]}, constrained_layout=True)
   axs = axs.flatten()
   f.suptitle('Time = {:.2f}'.format(current_time))
   
   #################################################
   # pressure image
   im_p = axs[0].contourf(pressure, extent=plenum_extent, **props(titles[0]))
   for tube_extent, tube_p in zip(tube_extents, Tube_p):
      midpoint = 0.5 * (tube_extent[2] + tube_extent[3])
      x = np.linspace(tube_extent[0], tube_extent[1], tube_n_cells)
      y_upper = midpoint + 0.5 * np.array([Area(xi) for xi in x]) * diameter_tube
      y_lower = midpoint - 0.5 * np.array([Area(xi) for xi in x]) * diameter_tube
      y_upper_scale = 0.2
      y_upper += y_upper_scale*np.max(y_upper) # increase y_upper by y_upper_scale (in percent) to avoid overlapping in the image
      lower = midpoint - 2.0 * diameter_tube
      upper = midpoint + 2.0 * diameter_tube
      tube_extent[2] = lower
      tube_extent[3] = upper
      im_p = axs[0].contourf(tube_p, extent=tube_extent, **props(titles[0]))
      axs[0].fill_between(x, y_upper, np.max(y_upper), color='white')
      axs[0].fill_between(x, y_lower, np.min(y_lower), color='white')
   axs[0].set_title(titles[0])
   cbar = plt.colorbar(im_p, ax=axs[0])
   cbar.set_label('[bar]', rotation=270, labelpad=15)
   
   #################################################
   # temperature image
   # im_T = axs[1].imshow(temperature * T_ref, extent=plenum_extent, **props(titles[1]))
   im_T = axs[1].contourf(temperature * T_ref, extent=plenum_extent, **props(titles[1]))
   axs[1].set_title(titles[1])
   # axs[1].set(aspect='equal')
   cbar = plt.colorbar(im_T, ax=axs[1])
   cbar.set_label('[K]', rotation=270, labelpad=15)
   
   
   ################################################
   # Passive Scalar
   
   im_X = axs[2].contourf(passiveScalarMF, extent=plenum_extent, **props(titles[2]))
   for tube_extent, tube_X in zip(tube_extents, Tube_X):
      midpoint = 0.5 * (tube_extent[2] + tube_extent[3])
      x = np.linspace(tube_extent[0], tube_extent[1], tube_n_cells)
      y_upper = midpoint + 0.5 * np.array([Area(xi) for xi in x]) * diameter_tube
      y_lower = midpoint - 0.5 * np.array([Area(xi) for xi in x]) * diameter_tube
      y_upper += y_upper_scale*np.max(y_upper) # increase y_upper by y_upper_scale (in percent) to avoid overlapping in the image
      im_X = axs[2].contourf(tube_X, extent=tube_extent, **props(titles[2]))
      axs[2].fill_between(x, y_upper, np.max(y_upper), color='white')
      axs[2].fill_between(x, y_lower, np.min(y_lower), color='white')
   axs[2].set_title(titles[2])
   cbar = plt.colorbar(im_X, ax=axs[2])

   
   #####################################################
   # velocity plot
   velU = rhoU / rho #* u_ref
   velV = rhoV / rho #* u_ref
   skip = 5
   # scale = 2.
   props_dict = props(titles[3])

   x = np.linspace(*plenum_extent[0:2], num=velU[::skip, ::skip].shape[1], endpoint=True)
   y = np.linspace(*plenum_extent[2:], num=velU[::skip, ::skip].shape[0], endpoint=True)
   passiveScalarMF,Y = np.meshgrid(x,y)
   # Q = axs[3].quiver(passiveScalarMF, Y, velU[::skip,::skip], velV[::skip,::skip], scale=u_ref, units='inches', width=0.01)
   Q = axs[3].quiver(passiveScalarMF, Y, velU[::skip,::skip], velV[::skip,::skip], **props_dict)
   axs[3].quiverkey(Q, 0.5, 1.05, 1./props_dict['scale'], '{}'.format(round(u_ref,1))+r'$ \frac{m}{s}$', labelpos='E', coordinates='axes')
   # axs[3].set_title('Velocity Field')
   # axs[3].set(aspect='equal')
   
   # f.subplots_adjust(hspace=0.4)
   f.savefig('{}/Figure{:05d}.png'.format(output_path, i), bbox_inches='tight', dpi=100)
   f.clear()
   plt.close(f)


# Call ffmpeg to make a movie

# os.system('ffmpeg -framerate 20 -i {}/Figure%05d.png -crf 20 {}/../Movie.mkv'.format(output_path, output_path))

# use start_number if not starting at File_00000.png
os.system('ffmpeg -start_number {} -framerate 20 -i {}/Figure%5d.png -crf 20 {}/../Movie.mkv'.format(t_index_array[0], output_path, output_path))