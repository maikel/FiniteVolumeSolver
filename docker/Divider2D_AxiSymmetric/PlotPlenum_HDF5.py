import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
sys.path.append(pathname)
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
   dataPath = 'ConvergentNozzle'


plenum = "{}/Plenum.h5".format(dataPath)
outPath = dataPath
output_path = '{}/Visualization/Plenum/'.format(outPath)

# Get nsteps from plenum
file = h5py.File(plenum, mode='r', swmr=True)
nsteps = file['data'].shape[0]
file.close()


def PrintProgress(i):
  progress = int(100.0 * float(i) / (nsteps - 1))
  print('[{:3d}%] Reading slice [{}/{}]'.format(progress, i + 1, nsteps))
os.makedirs(output_path, exist_ok=True)

import math

for i in range(nsteps):
   PrintProgress(i)
   
   plenum_variables = ["Pressure", "Density", "Momentum_0", "Momentum_1", "Temperature", "H2", 'vfrac']
   (p, rho, rhou, rhov, rhoH2, T, vols), current_time, extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   # # alternative:
   # plenum_data, current_time, extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   # p = plenum_data[plenum_dict['Pressure']]
   # rho = plenum_data[plenum_dict['Density']]
   # rhou = plenum_data[plenum_dict['Momentum_0']]
   # rhov = plenum_data[plenum_dict['Momentum_1']]
   # vols = plenum_data[plenum_dict['vfrac']]
   
   f, axs = plt.subplots(nrows=2, ncols=2, figsize=(35. / 2, 15. / 2), gridspec_kw={'width_ratios': [3,1]})
   axs = axs.flatten()
   f.suptitle('Time = {:.2f}'.format(current_time))
   # Pressure image
   p = np.ma.masked_array(p, vols < 1e-14)
   levels = np.linspace(0.9, 3.0, 30)
   pressure_options = {
     'origin': 'lower',
     'cmap': 'jet',
     'levels': levels, 
     'vmin': levels[0],
     'vmax': levels[-1],
     'extend': 'both'
   }
   im_p = axs[0].contourf(p / 101325.0, extent=extent, **pressure_options)
   axs[0].set_title('Pressure')
   axs[0].set(aspect='equal')
   cbar = plt.colorbar(im_p, ax=axs[0])
   cbar.set_label('[bar]', rotation=270, labelpad=15)
   
   # Temperature image
   T = np.ma.masked_array(T, vols < 1e-14)
   im_T = axs[1].imshow(T, origin='lower', vmin=290, vmax=3000.0, interpolation='none', extent=extent)
   axs[1].set_title('Temperature')
   axs[1].set(aspect='equal')
   cbar = plt.colorbar(im_T, ax=axs[1])
   cbar.set_label('[K]', rotation=270, labelpad=15)
   
   # Passive Scalar
   rho = np.ma.masked_array(rho, vols < 1e-14)
   H2 = rhoH2 / rho
   levels = np.linspace(0.0, 1.0, 40)
   im_X = axs[2].contourf(H2, origin='lower', extent=extent, levels=levels, extend='both', cmap='twilight')
   axs[2].set_title('H2 mass fractions')
   axs[2].set(aspect='equal')
   cbar = plt.colorbar(im_X, ax=axs[2])

   # velocity plot

   u_ref = math.sqrt(101325.)
   u = rhou / rho / u_ref
   v = rhov / rho / u_ref
   # print(np.max(u))
   skip = 5
   scale = 4.
   x = np.linspace(*extent[0:2], num=u[::skip, ::skip].shape[1], endpoint=True)
   y = np.linspace(*extent[2:], num=u[::skip, ::skip].shape[0], endpoint=True)
   X,Y = np.meshgrid(x,y)
   # Q = axs[3].quiver(X, Y, u[::skip,::skip], v[::skip,::skip], scale=u_ref, units='inches', width=0.01)
   Q = axs[3].quiver(X, Y, u[::skip,::skip], v[::skip,::skip], angles='xy', scale_units='xy', scale=scale)
   axs[3].quiverkey(Q, 0.5, 1.05, 1./scale, '{}'.format(round(u_ref,1))+r'$ \frac{m}{s}$', labelpos='E', coordinates='axes')
   # axs[3].set_title('Velocity Field')
   axs[3].set(aspect='equal')
   
   # f.subplots_adjust(hspace=0.4)
   f.savefig('{}/Figure{:04d}.png'.format(output_path, i), bbox_inches='tight')
   f.clear()
   plt.close(f)

os.system('ffmpeg -y -framerate 20 -i {}/Figure%04d.png -crf 20 {}/../Movie.mkv'.format(output_path, output_path))