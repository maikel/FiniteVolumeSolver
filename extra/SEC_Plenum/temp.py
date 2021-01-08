import yt
import os
import numpy as np
import matplotlib.pyplot as plt

plot_string = "quiver"
dataPath = "/home/christian/FVS_develop/FiniteVolumeSolver/build_Release_2d_EB/SEC_Plenum/"
plenum_number = 0
MASK_DATA_WITH_VFRAC = True
vfrac_cutoff = 1.0e-16

def loadQuiverData(plotfile):
  ds = yt.load(plotfile)
  ad = ds.all_data()
  current_time = ds.current_time
  rho = np.array(ad.index.grids[0]['Density']) #assume only one grid! --> TODO: Extend to possible more grids
  mom0 = np.array(ad.index.grids[0]['Momentum_0'])
  mom1 = np.array(ad.index.grids[0]['Momentum_1'])
  mom0 = np.squeeze( mom0 )
  mom1 = np.squeeze( mom1 )
  rho = np.squeeze( rho )
  if MASK_DATA_WITH_VFRAC:
    vfrac = np.array(ad.index.grids[0]["vfrac"])
    vfrac = np.squeeze(vfrac)
    rho = np.ma.masked_where( vfrac<=vfrac_cutoff, rho)
    mom0 = np.ma.masked_where( vfrac<=vfrac_cutoff, mom0)
    mom1 = np.ma.masked_where( vfrac<=vfrac_cutoff, mom1)
  return mom0, mom1, rho, current_time

def PrintProgress(i, plotfiles):
  ny = len(plotfiles)
  progress = int(100.0 * float(i) / (ny - 1))
  print('[{:3d}%] Reading plotfile {}'.format(progress, plotfiles[i]))

yt.funcs.mylog.setLevel(50)
output_dir = dataPath+'/Plotfiles/Plenum%i'%(plenum_number)

dirs = os.listdir(output_dir)
plotfiles = ['{}/{}'.format(output_dir, plt) for plt in dirs]
list.sort(plotfiles)
plotfiles = plotfiles[:]
# print(plotfiles)

ds = yt.load(plotfiles[0])
ad = ds.all_data()
LeftEgdes = ad.index.grids[0].LeftEdge
RightEdges = ad.index.grids[-1].RightEdge

m0, m1, rho, current_time = loadQuiverData(plotfiles[0])
PrintProgress(0, plotfiles)

nx = m0.shape
ny = len(plotfiles) # how many plotfiles do we have


ts = np.zeros((ny,1))
ts[0] = ds.current_time

if MASK_DATA_WITH_VFRAC:
  mom0 = np.ma.zeros(shape=(*nx, ny))
  mom1 = np.ma.zeros(shape=(*nx, ny))
  density = np.ma.zeros(shape=(*nx, ny))
else:
  mom0 = np.zeros(shape=(*nx, ny))
  mom1 = np.zeros_like(mom0)
  density = np.zeros_like(mom0)

mom0[:,:,0] = m0
mom1[:,:,0] = m1
density[:,:,0] = rho

output_path = dataPath+"Plotfiles/Plenum%i_%s/"%(plenum_number, plot_string)
if not os.path.exists(output_path):
  os.makedirs(output_path)

for i in range(1,ny):
   PrintProgress(i, plotfiles)
   m0, m1, rho, current_time = loadQuiverData(plotfiles[i])
   ts[i] = current_time
   mom0[:,:,i] = m0
   mom1[:,:,i] = m1
   density[:,:,i] = rho

def plotQuiver(plotfile, plot_string, mom0, mom1, density, current_time):
  fig, ax = plt.subplots(figsize=(10,10))
  # TODO: calculate x and y correct... something strange happen when we give ths values in ax.quiver...
  # x = np.linspace(LeftEgdes[0], RightEdges[0], num=mom0.shape[0], endpoint=True)
  # y = np.linspace(LeftEgdes[1], RightEdges[1], num=mom1.shape[1], endpoint=True)
  # X,Y = np.meshgrid(x,y)

  ux = mom0 / density
  uy = mom1 / density
  
  skip = 5
  scale = 1.0
  Q = ax.quiver(ux[::skip, ::skip], uy[::skip, ::skip], mom0[::skip, ::skip], 
        scale = scale, units='x', pivot='mid')
  # Q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], ux[::skip, ::skip], uy[::skip, ::skip], mom0[::skip, ::skip], 
  #      scale = scale, units='x', pivot='mid')

  cbar = fig.colorbar(Q)
  cbar.set_label("Momenttum0")
  ax.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
  ax.set(title="time = %2.6f"%current_time)


  fig.savefig(output_path+plot_string+"_"+plotfile.split("/")[-1]+".png", bbox='tight' )
  fig.clf()
  plt.close(fig)

for i in range(0, ny):
   plotQuiver(plotfiles[i], plot_string, mom0[:,:,i].T, mom1[:,:,i].T, density[:,:,i].T, ts[i])