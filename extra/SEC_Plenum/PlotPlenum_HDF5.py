import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
# print(FVS_path) 
sys.path.append(FVS_path+'/extra/')

import yt
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import h5py
import itertools


# dataPath = "/srv/public/Maikel/FiniteVolumeSolver/build_2D-Release/average_massflow"
# inputFilePath = "/srv/public/Maikel/FiniteVolumeSolver/examples/AMReX/EB/2D/"

dataPath = FVS_path+"/build_2D-Release/average_massflow"
inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"


# import importlib
# inputFile = importlib.import_module('SEC_Plenum')
# print(inputFile.y0s)
# print(inputFile.Area)
sys.path.append(inputFilePath)
from SEC_Plenum_Arrhenius import y0s, Area, tube_n_cells, p_ref, T_ref, rho_ref, Output
# print(y0s)
# print(Area)

plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization'.format(outPath)

def h5_load(path, num, variables):
   file = h5py.File(path, mode='r')
   strings = list(file['fields'].asstr())
   indices = [strings.index(var) for var in variables]
   shape = file['data'].shape
   nx = np.array([shape[2], shape[3]])
   datas = [np.reshape(np.array(file['data'][num, var, :, :]), [nx[1], nx[0]])  for var in indices]
   time = float(file['times'][num])
   dx = np.array(file['data'].attrs['cell_size'])
   xlower = np.array(file['data'].attrs['xlower'])
   xupper = xlower + dx * nx
   extent = [xlower[0], xupper[0], xlower[1], xupper[1]]
   file.close()
   return tuple(datas), time, extent


# Get nsteps from plenum
file = h5py.File(plenum, mode='r')
nsteps = file['data'].shape[0]
file.close()

# get data from output dictionary
output_dict = Output['outputs']
output_dict = [ out_dict for out_dict in output_dict if 'which_block' in out_dict] # strip 'CounterOutput'
for out_dict in output_dict:
   if 'Plenum' in out_dict['path']:
      if 'frequencies' in out_dict:
         plenum_out = int(out_dict['frequencies'][0])
      elif 'intervals' in out_dict:
         plenum_out = float(out_dict['intervals'][0])
      else:
         raise Exception("could not find output frequency or interval from plenum")
   elif 'Tube' in out_dict['path']: # assume all tubes have same frequency or intervals
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


def PrintProgress(i):
  progress = int(100.0 * float(i) / (nsteps - 1))
  print('[{:3d}%] Reading slice [{}/{}]'.format(progress, i + 1, nsteps))

os.makedirs(output_path, exist_ok=True)

tube_paths = ['{}/Tube{}.h5'.format(dataPath, tube) for tube in [0, 1 ,2]]

#for i, plotfile in itertools.dropwhile(lambda x: x[0] < 329, enumerate(plotfiles)):
#for i, plotfile in itertools.takewhile(lambda x: x[0] < 4, enumerate(plotfiles)):
for i in range(nsteps):
   PrintProgress(i)
   
   Tube_p = []
   extent_tubes = []
   
   for tube in tube_paths:
      (p_tube), current_time, extent = h5_load(tube, tube_output_factor*i, ["Pressure"])
      Tube_p.append(p_tube)
      extent_tubes.append(extent)
   
   def stackTubeDataTo2D(Tube_datalist):
     # all Tubedata is 1D but for contourf we need at least 2D data. so simply stack twice the 1d array
     for i, el in enumerate(Tube_datalist):
        el = np.squeeze(el)
        Tube_datalist[i] = np.stack((el,el))
     return Tube_datalist
   Tube_p = stackTubeDataTo2D(Tube_p)

   (p, rho, rhou, rhov, c, vols), current_time, extent = h5_load(plenum, i, ["Pressure", "Density", "Momentum_0", "Momentum_1", "SpeedOfSound", 'vfrac'])
   f, axs = plt.subplots(nrows=1, ncols=3, figsize=(40. / 2, 10. / 2), gridspec_kw={'width_ratios': [3, 1, 1]})
   f.suptitle('Time = {:.2f}'.format(current_time))
   # pressure image
   p = np.where(vols > 1e-14, p, np.nan)
   levels = np.linspace(1.9 * p_ref, 2.4 * p_ref, 30)
   pressure_options = {
     'origin': 'lower',
     'cmap': 'jet',
     'levels': levels,
     'vmin': levels[0],
     'vmax': levels[-1],
     'extend': 'both'
   }
   im_p = axs[0].contourf(p * p_ref, extent=extent, **pressure_options)
   for extent_tube, tube_p in zip(extent_tubes, Tube_p):
      midpoint = 0.5 * (extent_tube[2] + extent_tube[3])
      D = 0.03
      x = np.linspace(extent_tube[0], extent_tube[1], tube_n_cells)
      y_upper = midpoint + 0.5 * np.array([Area(xi) for xi in x]) * D
      y_lower = midpoint - 0.5 * np.array([Area(xi) for xi in x]) * D
      lower = midpoint - 4.0 * D
      upper = midpoint + 4.0 * D
      extent_tube[2] = lower
      extent_tube[3] = upper
      im_p = axs[0].contourf(tube_p * p_ref, extent=extent_tube, **pressure_options)
      axs[0].fill_between(x, y_upper, np.max(y_upper), color='white')
      axs[0].fill_between(x, y_lower, np.min(y_lower), color='white')
   axs[0].set_title('Pressure')
   axs[0].set(aspect='equal')
   plt.colorbar(im_p, ax=axs[0])
   # temperature image
   T = p / rho
   im_T = axs[1].imshow(T * T_ref, origin='lower', vmin=1.0 * T_ref, vmax=15 * T_ref, interpolation='none', extent=extent)
   axs[1].set_title('Temperature')
   axs[1].set(aspect='equal')
   plt.colorbar(im_T, ax=axs[1])
   # velocity field   rho = np.ma.masked_array(rho, vols > 1e-14)
   rho = np.ma.masked_array(rho, vols < 1e-14)
   u = rhou / rho
   v = rhov / rho
   skip = 5
   #scale = 40.0
   x = np.linspace(*extent[0:2], num=u[::skip, ::skip].shape[1], endpoint=True)
   y = np.linspace(*extent[2:], num=u[::skip, ::skip].shape[0], endpoint=True)
   X,Y = np.meshgrid(x,y)
   Q = axs[2].quiver(X, Y, u[::skip,::skip], v[::skip,::skip], scale=4, units='inches', width=0.01)
   axs[2].quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
   axs[2].set_title('Velocity Field')
   axs[2].set(aspect='equal')
   f.savefig('{}/Figure{:04d}.png'.format(output_path, i), bbox_inches='tight')
   # f.tight_layout()
   f.subplots_adjust(left=0.3)
   f.clear()
   plt.close(f)
