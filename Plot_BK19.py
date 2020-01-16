import os
import yt
import yt.funcs as funcs
import matplotlib as mpl
import matplotlib.pyplot as plt
import shutil
import numpy as np
import matplotlib.style
from mpl_toolkits.axes_grid1 import AxesGrid
#%matplotlib inline

if (not funcs.get_interactivity()):
  yt.toggle_interactivity()

def Load(path):
  #p = '{}/'.format(path)
  p = path
  shutil.copy2('extra/yt/WarpXHeader', p)
  shutil.copy2('extra/yt/warpx_job_info', p)
  return yt.load(p)

def PlotVariables(base_dir, step_dir, timestep):

  path = os.path.join(base_dir, step_dir, 'plt{:09d}'.format(timestep))
  ds   = Load(path)

  variables = ['Density', 'PTdensity', 'Velocity_0', 'Velocity_1']
  p = yt.plot_2d(ds, variables, origin='native')
  p.set_buff_size(256) # to remove some interpolation errors

  fig = plt.figure()
  for i, var in enumerate(variables):
    ax = plt.subplot(2,3,i+1)
    im = ax.imshow(np.array(p.frb[var]), origin='lower',extent=np.array(p.bounds))
    im.set_cmap('Blue-Red')
    plt.colorbar(im, ax=ax, shrink=.75)
    ax.set_title(var)

  g = ds.index.grids[0]
  data = np.squeeze(np.array(g['pi']))

  extent = [g.LeftEdge[0]  - g.dds[0]/2,
            g.RightEdge[0] + g.dds[0]/2,
            g.LeftEdge[1]  - g.dds[1]/2,
            g.RightEdge[1] + g.dds[1]/2]

  data_nd = np.zeros(g.shape[:2]+1)
  data_nd[:-1,:-1] = data[:,:,0]
  data_nd[:-1, -1] = data[:,-1,1]
  data_nd[ -1,:-1] = data[-1,:,2]
  data_nd[ -1, -1] = data[-1,-1,3]

  ax = plt.subplot(2,3,5)
  im = ax.imshow(data_nd, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('pi')

  fig.suptitle("timestep=%d, step=%s" %(timestep, step_dir))
  plt.show()

timestep = 1
base_dir = 'build2d'

step_dirs = ['BK19_Pseudo_Incompressible-pre-step',
             'BK19_Pseudo_Incompressible-advect',
             'BK19_Pseudo_Incompressible-advect-backward',
             'BK19_Pseudo_Incompressible-advect-backward-forward',
             'BK19_Pseudo_Incompressible-advect-backward-forward-advect',
             'BK19_Pseudo_Incompressible-advect-backward-forward-advect-backward']


PlotVariables(base_dir, step_dirs[0], timestep)
PlotVariables(base_dir, step_dirs[1], timestep)
PlotVariables(base_dir, step_dirs[2], timestep)
PlotVariables(base_dir, step_dirs[3], timestep)
PlotVariables(base_dir, step_dirs[4], timestep)
PlotVariables(base_dir, step_dirs[5], timestep)

#ds0 = yt.load(paths[0])
#ds1 = yt.load(paths[1])
#grid0 = ds0.covering_grid(level=0, left_edge=(0.0, 0.0, 0.0), dims=ds.domain_dimensions)
#grid1 = ds1.covering_grid(level=0, left_edge=(0.0, 0.0, 0.0), dims=ds.domain_dimensions)
#u0 = np.array(grid0['Velocity_0'])
#u1 = np.array(grid1['Velocity_0'])
#du = u1 - u0
#mpl.pyplot.imshow(du[:,:,0].transpose(), interpolation='none')
#mpl.pyplot.colorbar()

