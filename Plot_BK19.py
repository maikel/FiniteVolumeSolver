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

def Load(base_dir, step_dir, timestep):
  path = os.path.join(base_dir, step_dir, 'plt{:09d}'.format(timestep))
  shutil.copy2('extra/yt/WarpXHeader', path)
  shutil.copy2('extra/yt/warpx_job_info', path)
  return yt.load(path)

def LoadPv(base_dir, step_dir, timestep):
  path = os.path.join(base_dir, step_dir, 'Pv_plt{:09d}'.format(timestep))
  shutil.copy2('extra/yt/WarpXHeader', path)
  shutil.copy2('extra/yt/warpx_job_info', path)
  return yt.load(path)

def PlotVariables(base_dir, step_dir, timestep, fignum = -1):

  ds = Load(base_dir, step_dir, timestep)

  variables = ['Density', 'PTdensity', 'Velocity_0', 'Velocity_1']
  p = yt.plot_2d(ds, variables, origin='native')
  p.set_buff_size(256) # to remove some interpolation errors

  if (fignum == -1):
    fig = plt.figure()
  else:
    fig = plt.figure(fignum)
    fig.clf()

  for i, var in enumerate(variables):
    ax = plt.subplot(2,4,i+1)
    im = ax.imshow(np.array(p.frb[var]), origin='lower', extent=np.array(p.bounds))
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

  ax = plt.subplot(2,4,5)
  im = ax.imshow(data_nd.T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('pi')

  #ds = LoadPv(base_dir, step_dir, timestep)

  ##variables = ['Pv_0', 'Pv_1']
  #variables = ['Pv_0']
  #p = yt.plot_2d(ds, variables, origin='native')
  #p.set_buff_size(256) # to remove some interpolation errors

  #for i, var in enumerate(variables):
    #ax = plt.subplot(2,4,i+6)
    #im = ax.imshow(np.array(p.frb[var]), origin='lower', extent=np.array(p.bounds))
    #im.set_cmap('Blue-Red')
    #plt.colorbar(im, ax=ax, shrink=.75)
    #ax.set_title(var)

  ds = Load(base_dir, step_dir, timestep)
  g = ds.index.grids[0]

  data = np.squeeze(np.array(g['Pv_x']))

  extent = [g.LeftEdge[0]  - g.dds[0]/2,
            g.RightEdge[0] + g.dds[0]/2,
            g.LeftEdge[1],
            g.RightEdge[1]]

  data_fc = np.zeros((g.shape[0]+1, g.shape[1]))
  data_fc[:-1,:] = data[:,:,0]
  data_fc[-1,:] = data[-1,:,1]
  div = data[:,:,1] - data[:,:,0]

  ax = plt.subplot(2,4,6)
  im = ax.imshow(data_fc.T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('Pv_x')

  data = np.squeeze(np.array(g['Pv_y']))

  extent = [g.LeftEdge[0],
            g.RightEdge[0],
            g.LeftEdge[1]  - g.dds[1]/2,
            g.RightEdge[1] + g.dds[1]/2]

  data_fc = np.zeros((g.shape[0], g.shape[1]+1))
  data_fc[:,:-1] = data[:,:,0]
  data_fc[:,-1] = data[:,-1,1]
  div = div + data[:,:,1] - data[:,:,0]

  ax = plt.subplot(2,4,7)
  im = ax.imshow(data_fc.T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('Pv_y')

  ax = plt.subplot(2,4,8)
  im = ax.imshow(div.T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('div(Pv)')

  fig.suptitle("timestep=%d, t=%.4f, stage=%s" %(timestep, ds.current_time.value, step_dir[5:]))
  plt.show()

timestep = 1
base_dir = 'build2d'

step_dirs = ['BK19_pre-step',
             'BK19_advect',
             'BK19_advect-backward',
             'BK19_advect-backward-forward',
             'BK19_advect-backward-forward-advect',
             'BK19_advect-backward-forward-advect-backward']


PlotVariables(base_dir, step_dirs[0], timestep, 1)
PlotVariables(base_dir, step_dirs[1], timestep, 2)
PlotVariables(base_dir, step_dirs[2], timestep, 3)
PlotVariables(base_dir, step_dirs[3], timestep, 4)
PlotVariables(base_dir, step_dirs[4], timestep, 5)
PlotVariables(base_dir, step_dirs[5], timestep, 6)

#ds0 = yt.load(paths[0])
#ds1 = yt.load(paths[1])
#grid0 = ds0.covering_grid(level=0, left_edge=(0.0, 0.0, 0.0), dims=ds.domain_dimensions)
#grid1 = ds1.covering_grid(level=0, left_edge=(0.0, 0.0, 0.0), dims=ds.domain_dimensions)
#u0 = np.array(grid0['Velocity_0'])
#u1 = np.array(grid1['Velocity_0'])
#du = u1 - u0
#mpl.pyplot.imshow(du[:,:,0].transpose(), interpolation='none')
#mpl.pyplot.colorbar()

