workdir_path = '/srv/public/Maikel/FiniteVolumeSolver/extra/'

import sys
sys.path.append(workdir_path)

import yt
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import amrex.plotfiles

base_data_path = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Release/Divider_DE5_800'
plotfiles_path = '{}/Plotfiles'.format(base_data_path)
dirs = os.listdir(plotfiles_path)
plotfiles = ['{}/{}'.format(plotfiles_path, plt) for plt in dirs]
list.sort(plotfiles)
#plotfiles = [plotfiles[-1]]

figure_out_path = '{}/Divider/Visualization'.format(workdir_path)
os.makedirs(figure_out_path, exist_ok=True)

def CreateH5(path, plotfiles):
   import h5py
   vars = ['Density', 'Momentum_0', 'Momentum_1', 'Energy', 'Pressure', 'SpeedOfSound', 'vfrac']
   fields = [n.encode("ascii", "ignore") for n in vars]
   datas, current_time, extents = amrex.plotfiles.yt_load(plotfiles[0], vars)
   file = h5py.File(path, 'w')
   nsteps = len(plotfiles)
   nvars = len(vars)
   (nx, ny) = datas[0].shape
   data = file.create_dataset('/data', shape=(nvars, nx, ny, nsteps))
   dx = np.array([(extents[1] - extents[0]) / nx, (extents[3] - extents[2]) / ny])
   data.attrs['cell_size'] = dx
   data.attrs['xlower'] = extents[0]
   times = file.create_dataset('/times', shape=(nsteps))
   times[0] = current_time
   file.create_dataset('/fields', data=fields)
   for i, var in enumerate(datas):
      data[i, :, :, 0] = var
   for plotfile in plotfiles[1:]:
      datas, current_time, __ = amrex.plotfiles.yt_load(plotfiles[0], vars)
      times[0] = current_time
      for i, var in enumerate(datas):
         data[i, :, :, 0] = var
   file.close()

CreateH5('test.h5', plotfiles)

# for i, plotfile in enumerate(plotfiles):
#    (p, rho, vols), current_time, extent = amrex.plotfiles.yt_load(plotfile, ['Pressure', 'Density', 'vfrac'])
#    p = np.ma.masked_array(p, vols < 1e-14)
#    rho = np.ma.masked_array(rho, vols < 1e-14)
#    print('{}: {}, {}'.format(plotfile, current_time, extent))
#    f, axs = plt.subplots(nrows=2, ncols=1, figsize=(25, 20))
#    f.suptitle('Time = {:.2e}'.format(float(current_time)))
#    # pressure image
#    im_p = axs[0].imshow(p, origin='lower', vmin=0.9e5, vmax=np.max(p), interpolation='none', extent=extent, cmap='jet')
#    axs[0].set_title('Pressure')
#    plt.colorbar(im_p, ax=axs[0])
#    # temperature image
#    # im_rho = axs[1].imshow(rho, origin='lower', vmin=1.1, vmax=np.max(rho), interpolation='none', extent=extent, cmap='jet')
#    im_rho = axs[1].contourf(rho, origin='lower', levels=np.linspace(1., 3.0, 30), vmin=1., vmax=3.0, extent=extent, cmap='jet', extend='both')
#    axs[1].set_title('Density')
#    plt.colorbar(im_rho, ax=axs[1])
#    # save figure
#    plt.tight_layout()
#    figure_path = '{}/Figure{:04d}.png'.format(figure_out_path, i)
#    print('Save figure to {}.'.format(figure_path))
#    f.savefig(figure_path)
#    f.clear()
#    plt.close(f)