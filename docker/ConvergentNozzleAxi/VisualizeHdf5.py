import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
sys.path.append(pathname)
import amrex.h5_io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import h5py
import numpy as np

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

case = 'De1'
Ma = 1.61
nx = 1200
path = '{}/Divider_{}_nx-{}_Ma-{}.h5'.format(pathname, case, nx, Ma)
output_path = '{}/Figures/Hdf5/nx_{}_Mach_{}/'.format(pathname, nx, Ma)
os.makedirs(output_path, exist_ok=True)

vmin = None #0.4
vmax = None #2.0

file = h5py.File(path, mode='r')
nsteps = file['data'].shape[0]
file.close()


for i in range(nsteps):
    (rho, p, vfrac), time, extent, __ = amrex.h5_io.h5_load_spec_timepoint_variable(path, i, ['Density', 'Pressure', 'vfrac'])
    p = np.ma.masked_array(p, vfrac < 1e-14)
    rho = np.ma.masked_array(rho, vfrac < 1e-14)

    f, axs = plt.subplots(nrows=2, ncols=1, figsize=(25, 20))
    f.suptitle('Time = {:.2e}'.format(time))
    im_p = axs[0].imshow(p, origin='lower', interpolation='none', extent=extent, cmap='jet')
    axs[0].set_title('Pressure')
    plt.colorbar(im_p, ax=axs[0])
    # Temperature image
    Rspec = 287.058
    T = p / rho / Rspec
    im_T = axs[1].imshow(T, origin='lower', interpolation='none', extent=extent, cmap='jet')
    axs[1].set_title('Temperature')
    plt.colorbar(im_T, ax=axs[1])
    # save figure
    figure_path = '{}/Figure{:04d}.png'.format(output_path, i)
    print('Save figure to {}.'.format(figure_path))
    f.savefig(figure_path, bbox_inches='tight')
    f.clear()
    plt.close(f)
