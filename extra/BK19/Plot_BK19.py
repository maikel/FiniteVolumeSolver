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

#if (not funcs.get_interactivity()):
  #yt.toggle_interactivity()

def _vel0(field, data):
    return data['Momentum_0'] / data['Density']

def _vel1(field, data):
    return data['Momentum_1'] / data['Density']

def Load(base_dir, step_dir, timestep):
  path = os.path.join(base_dir, step_dir, 'partition_0_plt{:09d}'.format(timestep))
  shutil.copy2('../yt/WarpXHeader', path)
  shutil.copy2('../yt/warpx_job_info', path)

  ds = yt.load(path)
  ds.add_field(('boxlib', 'Velocity_0'), function=_vel0, units='', sampling_type='cell')
  ds.add_field(('boxlib', 'Velocity_1'), function=_vel1, units='', sampling_type='cell')

  return ds

def PlotVariables(ds, variables, label, fignum = -1):

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

  ax = plt.subplot(2,4,6)
  im = ax.imshow(data_nd.T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('pi')

  #ds = Load(base_dir, step_dir, timestep)
  #g = ds.index.grids[0]

  #data = np.squeeze(np.array(g['Pv_x']))

  #extent = [g.LeftEdge[0]  - g.dds[0]/2,
            #g.RightEdge[0] + g.dds[0]/2,
            #g.LeftEdge[1],
            #g.RightEdge[1]]

  #data_fc = np.zeros((g.shape[0]+1, g.shape[1]))
  #data_fc[:-1,:] = data[:,:,0]
  #data_fc[-1,:] = data[-1,:,1]
  #div = data[:,:,1] - data[:,:,0]

  #ax = plt.subplot(2,4,6)
  #im = ax.imshow(data_fc.T, origin='lower', extent=extent)
  #im.set_cmap('Blue-Red')
  #plt.colorbar(im, ax=ax, shrink=.75)
  #ax.set_title('Pv_x')

  #data = np.squeeze(np.array(g['Pv_y']))

  #extent = [g.LeftEdge[0],
            #g.RightEdge[0],
            #g.LeftEdge[1]  - g.dds[1]/2,
            #g.RightEdge[1] + g.dds[1]/2]

  #data_fc = np.zeros((g.shape[0], g.shape[1]+1))
  #data_fc[:,:-1] = data[:,:,0]
  #data_fc[:,-1] = data[:,-1,1]
  #div = div + data[:,:,1] - data[:,:,0]

  #ax = plt.subplot(2,4,7)
  #im = ax.imshow(data_fc.T, origin='lower', extent=extent)
  #im.set_cmap('Blue-Red')
  #plt.colorbar(im, ax=ax, shrink=.75)
  #ax.set_title('Pv_y')

  #ax = plt.subplot(2,4,8)
  #im = ax.imshow(div.T, origin='lower', extent=extent)
  #im.set_cmap('Blue-Red')
  #plt.colorbar(im, ax=ax, shrink=.75)
  #ax.set_title('div(Pv)')

  fig.suptitle("timestep=%d, t=%.4f, stage=%s" %(timestep, ds.current_time.value, label))
  plt.show()

timestep = 1
base_dir = '../../build2d/Debug'

substeps = ['BK19_pre-step',
            'BK19_advect',
            'BK19_advect-backward',
            'BK19_advect-backward-forward',
            'BK19_advect-backward-forward-advect',
            'BK19_advect-backward-forward-advect-backward']

vars = ['Density', 'Velocity_0', 'Velocity_1', 'PTdensity', 'PTinverse']
#vars = ['Density', 'Momentum_0', 'Momentum_1', 'PTdensity', 'PTinverse']

for i in range(len(substeps)):
  ds = Load(base_dir, substeps[i], timestep)
  ds.field_info[('boxlib', 'Momentum_0')].take_log = False
  ds.field_info[('boxlib', 'Momentum_1')].take_log = False
  ds.field_info[('boxlib', 'Velocity_0')].take_log = False
  ds.field_info[('boxlib', 'Velocity_1')].take_log = False
  PlotVariables(ds, vars, substeps[i][5:], i+1)
  #PlotVariables(ds, vars, 'step = {0:3}'.format(timestep), i+1)

ds = Load(base_dir, 'BK19_advect-backward', timestep)

g = ds.index.grids[0]

fig = plt.figure(10)
fig.clf()

data = np.squeeze(np.array(g['Pu_faces']))

extent = [g.LeftEdge[0]  - g.dds[0]/2,
          g.RightEdge[0] + g.dds[0]/2,
          g.LeftEdge[1],
          g.RightEdge[1]]

data_fc = np.zeros((g.shape[0]+1, g.shape[1]))
data_fc[:-1,:] = data[:,:,0]
data_fc[-1,:] = data[-1,:,1]

ax = plt.subplot(1,2,1)
im = ax.imshow(data_fc.T, origin='lower', extent=extent)
im.set_cmap('Blue-Red')
plt.colorbar(im, ax=ax, shrink=.75)
ax.set_title('Pu_faces')

data = np.squeeze(np.array(g['Pv_faces']))

extent = [g.LeftEdge[0],
          g.RightEdge[0],
          g.LeftEdge[1]  - g.dds[1]/2,
          g.RightEdge[1] + g.dds[1]/2]

data_fc = np.zeros((g.shape[0], g.shape[1]+1))
data_fc[:,:-1] = data[:,:,0]
data_fc[:,-1] = data[:,-1,1]

ax = plt.subplot(1,2,2)
im = ax.imshow(data_fc.T, origin='lower', extent=extent)
im.set_cmap('Blue-Red')
plt.colorbar(im, ax=ax, shrink=.75)
ax.set_title('Pv_faces')

plt.show()
