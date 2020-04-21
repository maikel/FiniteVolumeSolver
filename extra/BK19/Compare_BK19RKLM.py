import os
import h5py
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

def LoadFVS(base_dir, step_dir, partname, timestep):
  path = os.path.join(base_dir, step_dir, partname+'plt{:09d}'.format(timestep))
  shutil.copy2('../yt/WarpXHeader', path)
  shutil.copy2('../yt/warpx_job_info', path)
  return yt.load(path)

def LoadRKLM(base_dir, filename, ts_label):
  inner = (slice(2,-2), slice(2,-2))
  path = os.path.join(base_dir, filename)
  file = h5py.File(path,'r')

  data = {}
  for var in file.keys():
    tmp = file[var].get(var+'_'+ts_label)
    if (tmp != None):
      data[var] = np.array(tmp[()])[inner]
  file.close()

  data['velu'] = data['rhou'] / data['rho']
  data['velv'] = data['rhov'] / data['rho']
  data['chi']  = 1.0 / data['Y']

  return data


def PlotVarComp(dsFVS, dsRKLM, vars_FVS, vars_RKLM, label, fignum = -1):

  p = yt.plot_2d(dsFVS, vars_FVS, origin='native')
  p.set_buff_size(64) # to remove some interpolation errors

  if (fignum == -1):
    fig = plt.figure()
  else:
    fig = plt.figure(fignum)
    fig.clf()

  numvars = len(vars_FVS)
  for i in range(numvars):
    ax = plt.subplot(3,numvars+1,1+i)
    im = ax.imshow(np.array(p.frb[vars_FVS[i]]), origin='lower', extent=np.array(p.bounds))
    im.set_cmap('Blue-Red')
    plt.colorbar(im, ax=ax, shrink=.75)
    ax.set_title(vars_FVS[i])

    ax = plt.subplot(3,numvars+1,(numvars+1)+1+i)
    im = ax.imshow(dsRKLM[vars_RKLM[i]].T, origin='lower', extent=np.array(p.bounds))
    im.set_cmap('Blue-Red')
    plt.colorbar(im, ax=ax, shrink=.75)
    ax.set_title(vars_RKLM[i])

    ax = plt.subplot(3,numvars+1,(2*(numvars+1))+1+i)
    im = ax.imshow(np.array(p.frb[vars_FVS[i]])-dsRKLM[vars_RKLM[i]].T, origin='lower', extent=np.array(p.bounds))
    im.set_cmap('Blue-Red')
    plt.colorbar(im, ax=ax, shrink=.75)
    ax.set_title(vars_FVS[i]+' - '+vars_RKLM[i])

  g = dsFVS.index.grids[0]
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

  ax = plt.subplot(3,numvars+1,numvars+1)
  im = ax.imshow(data_nd.T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('pi')

  ax = plt.subplot(3,numvars+1,(numvars+1)+numvars+1)
  im = ax.imshow(dsRKLM['p2_nodes'].T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('pi')

  ax = plt.subplot(3,numvars+1,2*(numvars+1)+numvars+1)
  im = ax.imshow(data_nd.T-dsRKLM['p2_nodes'].T, origin='lower', extent=extent)
  im.set_cmap('Blue-Red')
  plt.colorbar(im, ax=ax, shrink=.75)
  ax.set_title('difference')

  fig.suptitle(label)
  plt.show()

timestep = 1

basedir_FVS  = '../../build/Debug'
#basedir_FVS  = '../../build2d'
substeps_FVS = ['BK19_pre-step',
                'BK19_advect',
                'BK19_advect-backward',
                'BK19_advect-backward-forward',
                'BK19_advect-backward-forward-advect',
                'BK19_advect-backward-forward-advect-backward']
vars_FVS     = ['Density', 'Velocity_0', 'Velocity_1', 'PTdensity', 'PTinverse']
#vars_FVS     = ['Density', 'Momentum_0', 'Momentum_1', 'PTdensity', 'PTinverse']
partname = 'partition_0_'
#partname = ''

basedir_RKLM  = '/home/ray/git-projects/RKLM_Reference/output_travelling_vortex/'
filename_RKLM = "output_travelling_vortex_ensemble=1_64_64_1.0_comp_noip.h5"
substeps_RKLM = ['before_advect',
                 'after_advect',
                 'after_ebnaimp',
                 'after_efna',
                 'after_full_advect',
                 'after_full_ebnaexp']
vars_RKLM = ['rho', 'velu', 'velv', 'rhoY', 'chi']
#vars_RKLM = ['rho', 'rhou', 'rhov', 'rhoY', 'chi']

for i in range(len(substeps_FVS)):
  dsFVS  = LoadFVS(basedir_FVS, substeps_FVS[i], partname, timestep)
  dsRKLM = LoadRKLM(basedir_RKLM, filename_RKLM, '{0:03}_{1:s}'.format(timestep-1, substeps_RKLM[i]))
  PlotVarComp(dsFVS, dsRKLM, vars_FVS, vars_RKLM, substeps_FVS[i][5:]+', '+substeps_RKLM[i], i+1)
