import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
import amrex.h5_data_processing as dataManip
from amrex.other import import_file_as_module

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import h5py
import itertools

os.environ['HDF5_USE_FILE_LOCKING'] = 'True'

# check cli
if len(sys.argv)<2:
  errMsg = ('Not enough input arguments!\n'
               +'\t1. argument must be dataPath!\n'
               +'\toptional argument is name of the inputfile: --config=inputfile.py\n'
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

import_file_as_module(os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import Area, tube_n_cells, p_ref, rho_ref, Output, u_ref, t_ref
from inputfile import D as diameter_tube

try:
  from inputfile import T_ref, n_tubes  
except:
  from inputfile import T_ref
  n_tubes=1

plenum = "{}/Plenum.h5".format(dataPath)
outPath = '{}/Visualization/'.format(dataPath)
output_path = '{}/PlenumStartup/'.format(outPath)

# Get times and nSteps from plenum
print('Read times from {}'.format(plenum))
times = io.h5_load_timepoints(plenum)
nSteps = times.shape[0]
print('Found {} timepoints/steps'.format(nSteps))

# optional slicing in time-dimension
tplotmin = 0.
tplotmax = 160.
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

print('LOad passive scalar limits for scaling right the contourf plot.')
scalar_limits = io.getPassiveScalarLimits(plenum, t_index_array[0], t_index_array[-1])

titles = ['Pressure', 'Temperature',  'Passive Scalar', 'Velocity']
def props(title):
   nlevels = 101
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
      minT = 7.*T_ref
      maxT = 12*T_ref
      scaleT = 25
      nlevelsT = 61 #int(maxT-minT)//scaleT
      levels = np.linspace(minT, maxT, nlevelsT )
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
      levels = np.linspace(2.0, 6.0, nlevels)
      props.update(
         {
         'levels': levels,
         'vmin': levels[0],
         'vmax': levels[-1],
         'cmap': 'jet'
         }
      )
   if title == 'Species':
      levels = np.linspace(0., 1., nlevels)
      props.update(
         {
         'levels': levels,
         'vmin': levels[0],
         'vmax': levels[-1],
         'cmap': 'gray_r',
         'extend': 'neither'
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

# dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Pressure', 1.0, 8, output_path, True, PARALLEL)
# dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Density', 1.0,  8, output_path, True, PARALLEL)
# dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'Temperature', 1.0,  8, output_path, True, PARALLEL)
# dataManip.printSimpleStatsPlenumSingleTimepoint(np.zeros(2), 'PassiveScalar', 1.0,  8, output_path, True, PARALLEL)

def plotFigure(i, PARALLEL=False):
   if PARALLEL:
      print(f'plotting picture {i}')
   else:
      PrintProgress(i)
   Tube_p = []
   Tube_X = []
   Tube_fuel = []
   Tube_temperature = []
   tube_extents = []
   tube_variables = ["Pressure", "PassiveScalars", "Density", "Species"]

   for tube in tube_paths:
      tube_data, current_time, ext, tube_dict = io.h5_load_spec_timepoint_variable(tube, tube_output_factor*i, tube_variables)
      Tube_p.append(tube_data[tube_dict['Pressure']])
      Tube_X.append( tube_data[tube_dict['PassiveScalars']] / tube_data[tube_dict['Density']] )
      Tube_fuel.append( tube_data[tube_dict['Species']] / tube_data[tube_dict['Density']] )
      Tube_temperature.append(tube_data[tube_dict['Pressure']] / tube_data[tube_dict['Density']])
      tube_extents.append(ext)
   
   Tube_p = stackTubeDataTo2D(Tube_p)
   Tube_X = stackTubeDataTo2D(Tube_X)
   Tube_fuel = stackTubeDataTo2D(Tube_fuel)
   Tube_temperature = stackTubeDataTo2D(Tube_temperature)
   
   plenum_variables = ["Pressure", "Density", "Momentum_0", "Momentum_1", "PassiveScalars", 'vfrac']
   # # alternative:
   # (pressure, rho, rhoU, rhoV, rhoX, vols), current_time, plenum_extent, plenum_dict = io.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   plenum_data, current_time, plenum_extent, plenum_dict = io.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   volume_fraction = plenum_data[plenum_dict['vfrac']]

   volume_fraction_offset = 1.e-14
   pressure = np.ma.masked_array(plenum_data[plenum_dict['Pressure']], volume_fraction < volume_fraction_offset)
   rho = np.ma.masked_array(plenum_data[plenum_dict['Density']], volume_fraction < volume_fraction_offset)
   rhoU = np.ma.masked_array(plenum_data[plenum_dict['Momentum_0']], volume_fraction < volume_fraction_offset)
   rhoV = np.ma.masked_array(plenum_data[plenum_dict['Momentum_1']], volume_fraction < volume_fraction_offset)
   passiveScalarMF = np.ma.masked_array(plenum_data[plenum_dict['PassiveScalars']], volume_fraction < volume_fraction_offset) / rho
   temperature = pressure / rho

   # # # print out the first occurence of min/max value 
   # dataManip.printSimpleStatsPlenumSingleTimepoint(pressure, 'Pressure', current_time, output_path=output_path, PARALLEL=PARALLEL)
   # dataManip.printSimpleStatsPlenumSingleTimepoint(rho, 'Density', current_time, output_path=output_path, PARALLEL=PARALLEL)
   # dataManip.printSimpleStatsPlenumSingleTimepoint(temperature, 'Temperature', current_time, output_path=output_path, PARALLEL=PARALLEL)
   # dataManip.printSimpleStatsPlenumSingleTimepoint(passiveScalarMF, 'PassiveScalar', current_time, output_path=output_path, PARALLEL=PARALLEL)


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
      y_lower = midpoint - 0.75 * np.array([Area(xi) for xi in x]) * diameter_tube # nobody knows why 0.75 works better as 0.5 ...
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

   from mpl_toolkits.axes_grid1.inset_locator import inset_axes
   # inset axes....
   axins_fuel = axs[0].inset_axes([0.05, 0.6, 0.4, 0.3])
   im_fuel = axins_fuel.contourf(Tube_fuel[0], extent=tube_extent, **props('Species'))
   axins_fuel.fill_between(x, y_upper, np.max(y_upper), color='white')
   axins_fuel.fill_between(x, y_lower, np.min(y_lower), color='white')
   # sub region of the original image
   x1, x2, y1, y2 = tube_extent[0], -0.5, lower, upper
   axins_fuel.set_xlim(x1, x2)
   axins_fuel.set_ylim(y1, y2)
   axins_fuel.set_xticklabels([])
   axins_fuel.set_yticklabels([])
   axs[0].indicate_inset_zoom(axins_fuel, edgecolor="black")
   cax = inset_axes(axins_fuel,
                 width="5%",  # width = 10% of parent_bbox width
                 height="100%",  # height : 50%
                 loc='lower left',
                 bbox_to_anchor=(1.05, 0., 1, 1),
                 bbox_transform=axins_fuel.transAxes,
                 borderpad=0,
                 )
   # cax = axs[0].inset_axes([0.05+0.4, 0.5, 0.05, 0.4])
   cbar = plt.colorbar(im_fuel, cax=cax, ticks=[np.linspace(0.,1.,5)])
   cbar.set_label('Fuel [-]', rotation=270, labelpad=15)
  
   # plot in ax3 only upper half of the tube
   tube_extent_insetpressure = tube_extent.copy()
   tube_extent_insetpressure[2] = 0.

   # plot in ax4 the lower part
   tube_extent_insettemperature = tube_extent.copy()
   tube_extent_insettemperature[3] = 0.
   
   # inset axes....
   axins_pressure = axs[0].inset_axes([0.05, 0.1, 0.4, 0.3])
   axins_pressure.contourf(Tube_p[0], extent=tube_extent_insetpressure, **props('Pressure'))
   axins_pressure.contourf(Tube_temperature[0]*T_ref, extent=tube_extent_insettemperature, **props('Temperature'))
   
   axins_pressure.fill_between(x, y_upper, np.max(y_upper), color='white')
   axins_pressure.fill_between(x, y_lower, np.min(y_lower), color='white')
   # sub region of the original image
   x1, x2, y1, y2 = tube_extent[0], -0.5, lower, upper
   axins_pressure.set_xlim(x1, x2)
   axins_pressure.set_ylim(y1, y2)
   #axins_pressure.set_xticklabels([])
   axins_pressure.set_yticklabels([])
   axs[0].indicate_inset_zoom(axins_pressure, edgecolor="black")
   
   
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

if PARALLEL:
   values = ((i,PARALLEL) for i in t_index_array)

   with Pool(Num_CPUs) as pool:
      pool.starmap(plotFigure, values)
else:
   for i in t_index_array:
      plotFigure(i)
      exit()

# sort unorderd data and save them to disk
if PARALLEL:
   dataManip.sortSimpleStatsPlenum(outPath)

# Call ffmpeg to make a movie
try:
   # os.system('ffmpeg -framerate 20 -i {}/Figure%05d.png -crf 20 {}/../Movie.mkv'.format(output_path, output_path))

   # use start_number if not starting at File_00000.png
   os.system('ffmpeg -start_number {} -framerate 20 -i {}/Figure%5d.png -crf 20 {}/../PlenumMovieStartup.mkv'.format(t_index_array[0], output_path, output_path))
except:
   print('ffmpeg could not be started')
