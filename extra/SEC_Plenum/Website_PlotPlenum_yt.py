import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.yt_io as io
import amrex.other as other

import yt
import numpy as np
import itertools

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# check cli
if len(sys.argv)<2:
  errMsg = ('Not enough input arguments!\n'
               +'\t1. argument must be dataPath!\n'
               +'\toptional argument is name of the inputfile\n'
               +'\te.g. {} path --config=inputfile.py'.format(sys.argv[0])
            )
  raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

# name of the inputfile is optional
optional = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--config=' in el ]
if not optional:
    optional = ['inputfile.py'] # default value 
inputfileName = optional[0]

#'Divider2D_MassFlow_AxiSymmetric_Multiblock.py'

other.import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import n_level, d_tube


yt.funcs.mylog.setLevel(50)

output_dir = dataPath+'/Plotfiles/Plenum0'

dirs = os.listdir(output_dir)
plotfiles = ['{}/{}'.format(output_dir, plt) for plt in dirs]
list.sort(plotfiles)

output_path = '{}/Visualization'.format(dataPath)
os.makedirs(output_path, exist_ok=True)

def _velMag(field, data):
   velX = data["Momentum_0"]/data["Density"]
   velY = data["Momentum_1"]/data["Density"]

   magVel = np.sqrt( np.power(velX, 2) + np.power(velY, 2) )
   return magVel/330.

yt.add_field(
   name=('velocityMagnitude'),
   function=_velMag,
   take_log=False
)

plotfiles = [plotfiles[1005], plotfiles[1011]]
#plotfiles = [plotfiles[0]]

# f = plt.figure()
f, axs = plt.subplots(nrows=1, ncols=len(plotfiles), sharey=True, figsize=(4.5*len(plotfiles),2.5))
try:
   axs = axs.flatten()
except:
   axs = [axs]

text = ['a)',' b)']

# for i, plotfile in itertools.dropwhile(lambda x: x[0] < 1000, enumerate(plotfiles)):
#for i, plotfile in itertools.takewhile(lambda x: x[0] < 4, enumerate(plotfiles)):
# for i, plotfile in enumerate(plotfiles):
for i, plotfile in enumerate(plotfiles):

   ds = yt.load(plotfile)
   name = plotfile.split("/")[-1]
   print('Reading plotfile: {}'.format(name))

   (p, rho, rhou, rhov, c, vols), current_time, extent = io.yt_load(plotfile, ["Pressure", "Density", "Momentum_0", "Momentum_1", "SpeedOfSound", 'vfrac'])
   # possible values are: ['Density', 'Energy', 'Momentum_0', 'Momentum_1', 'Pressure', 'Species_0', 'Species_1', 'SpeedOfSound', 'vfrac']
   
   uVel = rhou / rho #np.where(vols > 1e-14, , np.nan)
   vVel = np.where(vols > 1e-14, rhov / rho, np.nan)
   magVel = np.sqrt(uVel**2+vVel**2)/330.

   # normalize extent with tube diameter
   extent /=d_tube
   # shrink extent
   extent_shrinked = [-0.5, 1.5, 0., 1.5 ]
   
   levels = np.linspace(0.0, 1.6, 16*2+1)
   velocity_options = {
   #   'origin': 'lower',
   #   'interpolation': 'hermite',
     'cmap': 'gist_heat_r',
     'vmin': levels[0],
     'vmax': levels[-1],
     'levels': levels,
     'extend': 'max',
   }

   if name.endswith('000'):
      extent_shrinked = [-0.5, 1.0, 0., 1.0 ]
      pass
   else:
      im_magvel = axs[i].contourf(magVel, extent=extent, **velocity_options)
      # axs[i].set_title('Time = {:.5f}'.format(current_time))
      if axs[i]==axs[-1]:
         cbar = f.colorbar(im_magvel, ax=axs.ravel().tolist(), 
                  shrink=1.0)
         cbar.set_label(r'$|\mathbf{u}|/a_0$')
      
      axs[i].annotate(text[i], xy=(0.05,0.85), xycoords='axes fraction', fontweight='bold', zorder=10)
   
   colors = ['k', 'b', 'g', 'r']
   colorIndex = 0
   # get cell size x from coarsest grid
   dds = ds.index.grids[0].dds[0] /d_tube

   # collectedGrids = io.yt_collect_child_grids(ds.index.grids[0])
   for grid in other.progressBar(ds.index.grids):
      # print(grid.dds)
      if dds > grid.dds[0]/d_tube:
         dds = grid.dds[0]/d_tube
         colorIndex +=1

      # draw level contour
      ledge = grid.LeftEdge /d_tube
      redge = grid.RightEdge /d_tube
      width = redge[0] - ledge[0]
      heigth = redge[1] - ledge[1]
      rect = patches.Rectangle(xy=ledge, width=width, height=heigth, fill=False, color=colors[colorIndex], lw=1.25, zorder=2)
      axs[i].add_patch(rect)
      
      # add rect Patches for grid with dds
      ddsX = grid.dds[0]/d_tube
      ddsY = grid.dds[1]/d_tube

      # print(grid.child_mask.shape)

      for nx in range(grid.child_mask.shape[0]):
         xedge = ledge[0] + nx * ddsX
         yedge = ledge[1]
         width = ddsX
         height = grid.child_mask.shape[1]*ddsY
         axs[i].add_patch(
            patches.Rectangle(
               xy=(xedge, yedge), 
               width=width, 
               height=height, 
               fill=False, 
               color='k',
               alpha = 0.15,
               lw=0.25,
               zorder=1
               )
         )
      for ny in range(grid.child_mask.shape[1]):
         xedge = ledge[0] 
         yedge = ledge[1] + ny * ddsY
         width = grid.child_mask.shape[0]*ddsX
         height = ddsY
         axs[i].add_patch(
            patches.Rectangle(
               xy=(xedge, yedge), 
               width=width, 
               height=height, 
               fill=False, 
               color='k',
               alpha = 0.15, 
               lw=0.25,
               zorder=1
               )
         )

   
   # white patch to hide the cutcells
   axs[i].add_patch(
      patches.Rectangle(
         xy=(extent[0], extent[2]+d_tube/d_tube/2.),
            width=(0.-extent[0]),
            height=(extent[3]-d_tube/d_tube/2.), 
            fill=True, color='w', zorder=2
         )
   )

   # patch at the cutcell boundary
   axs[i].add_patch(
      patches.Rectangle(
         xy=(extent[0], extent[2]+d_tube/d_tube/2.),
            width=(0.-extent[0]),
            height=(extent[3]-d_tube/d_tube/2.), 
            fill=False, color=colors[n_level-1], zorder=2
         )
   )

   axs[i].set(
      xlabel = r'$X/D$',
      ylabel = r'$Y/D$',
      aspect = 'equal',
      xlim = (extent_shrinked[:2]),
      ylim = (extent_shrinked[2:])
      )
   axs[i].label_outer()
   

f.savefig('{}/{}.png'.format(output_path, name), bbox_inches='tight', dpi=175)
f.savefig('{}/{}.pdf'.format(output_path, name), bbox_inches='tight', dpi=600)
f.clear()
plt.close(f)