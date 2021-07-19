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
import math

gamma = 1.4
Rspec = 287.058


os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1])
   inputFilePath = dataPath
else:
   sys.exit('ERROR: Please parse the dataPath from command Line!\n'
            +'try again as:\n'
            +'python3 {} dataPath'.format(sys.argv[0])
            )

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


for i in range(nsteps):
   PrintProgress(i)
   
   plenum_variables = ["Pressure", "Density", "Momentum_0", "Momentum_1", 'vfrac']
   (p, rho, rhou, rhov, vols), current_time, extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   
   p = np.ma.masked_array(p, vols < 1e-14)
   rho = np.ma.masked_array(rho, vols < 1e-14)
   rhou = np.ma.masked_array(rhou, vols < 1e-14)
   rhov = np.ma.masked_array(rhov, vols < 1e-14)
   u = rhou / rho
   v = rhov / rho

   # # alternative:
   # plenum_data, current_time, extent, plenum_dict = da.h5_load_spec_timepoint_variable(plenum, i, plenum_variables)
   # p = plenum_data[plenum_dict['Pressure']]
   # rho = plenum_data[plenum_dict['Density']]
   # rhou = plenum_data[plenum_dict['Momentum_0']]
   # rhov = plenum_data[plenum_dict['Momentum_1']]
   # vols = plenum_data[plenum_dict['vfrac']]
   
   f, axs = plt.subplots(nrows=2, ncols=2, figsize=(35. / 2, 15. / 2), gridspec_kw={'width_ratios': [3,3]})
   axs = axs.flatten()
   f.suptitle('Time = {:.2f}'.format(current_time))
   # Pressure image
   levels = np.linspace(0.9, 1.5, 40)
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
   
   # Density image
   im_rho = axs[1].imshow(rho, origin='lower', vmin=1.0, vmax=1.4, interpolation='none', extent=extent)
   axs[1].set_title('Density')
   axs[1].set(aspect='equal')
   cbar = plt.colorbar(im_rho, ax=axs[1])
   cbar.set_label('[kg/m^3]', rotation=270, labelpad=15)
   
   # Velocity X
   levels = np.linspace(0.0, 100., 40)
   im_u = axs[2].contourf(u, origin='lower', extent=extent, levels=levels, extend='both', cmap='twilight')
   axs[2].set_title('Velocity_X')
   axs[2].set(aspect='equal')
   cbar = plt.colorbar(im_u, ax=axs[2])
   cbar.set_label('[m/s]', rotation=270, labelpad=15)

   # Velocity Y
   levels = np.linspace(0.0, 100., 40)
   im_v = axs[3].contourf(v, origin='lower', extent=extent, levels=levels, extend='both', cmap='twilight')
   axs[3].set_title('Velocity_Y')
   axs[3].set(aspect='equal')
   cbar = plt.colorbar(im_v, ax=axs[3])
   cbar.set_label('[m/s]', rotation=270, labelpad=15)


   
   # f.subplots_adjust(hspace=0.4)
   f.savefig('{}/Figure{:04d}.png'.format(output_path, i), bbox_inches='tight')
   f.clear()
   plt.close(f)

os.system('ffmpeg -y -framerate 20 -i {}/Figure%04d.png -crf 20 {}/../Movie.mkv'.format(output_path, output_path))