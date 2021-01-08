workdir_path = '/srv/public/Maikel/FiniteVolumeSolver/extra/'

import sys
sys.path.append(workdir_path)

import yt
import os
import numpy as np
import matplotlib.pyplot as plt
import amrex.plotfiles

base_data_path = '/group/ag_klima/SFB1029_C01/Divider/Divider_DE5_2400'
plotfiles_path = '{}/Plotfiles'.format(base_data_path)
dirs = os.listdir(plotfiles_path)
plotfiles = ['{}/{}'.format(plotfiles_path, plt) for plt in dirs]
list.sort(plotfiles)
#plotfiles = [plotfiles[-1]]

figure_out_path = '{}/Divider/Visualization'.format(workdir_path)
os.makedirs(figure_out_path, exist_ok=True)

for i, plotfile in enumerate(plotfiles):
   (p, rho, vols), current_time, extent = amrex.plotfiles.load(plotfile, ['Pressure', 'Density', 'vfrac'])
   print('{}: {}, {}'.format(plotfile, current_time, extent))
   f, axs = plt.subplots(nrows=2, ncols=1, figsize=(25, 20))
   f.suptitle('Time = {:.2e}'.format(float(current_time)))
   # pressure image
   im_p = axs[0].imshow(p, origin='lower', vmin=0.9e5, vmax=1.4e5, interpolation='none', extent=extent, cmap='jet')
   axs[0].set_title('Pressure')
   plt.colorbar(im_p, ax=axs[0])
   # temperature image
   im_rho = axs[1].imshow(rho, origin='lower', vmin=1.1, vmax=1.5, interpolation='none', extent=extent, cmap='jet')
   axs[1].set_title('Density')
   plt.colorbar(im_rho, ax=axs[1])
   # save figure
   plt.tight_layout()
   figure_path = '{}/Figure{:04d}.png'.format(figure_out_path, i)
   print('Save figure to {}.'.format(figure_path))
   f.savefig(figure_path)
   f.clear()
   plt.close(f)