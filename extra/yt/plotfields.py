import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yt
import yt.funcs as funcs
import shutil

datapath = './BK19_Pseudo_Incompressible/plt000000000'

shutil.copy2('WarpXHeader', datapath)
shutil.copy2('warpx_job_info', datapath)

if (not funcs.get_interactivity()):
  yt.toggle_interactivity()

ds = yt.load( datapath )

ds.field_list

g = ds.index.grids[0]
pi = np.squeeze(np.array(g["pi"]))

x_nd = np.linspace(g.LeftEdge[0], g.RightEdge[0], g.shape[0]+1)
y_nd = np.linspace(g.LeftEdge[1], g.RightEdge[1], g.shape[1]+1)

a_nd = np.zeros(g.shape[:2]+1)
a_nd[:-1,:-1] = pi[:,:,0]
a_nd[:-1,-1] = pi[:,-1,1]
a_nd[-1,:-1] = pi[-1,:,2]
a_nd[-1,-1] = pi[-1,-1,3]

p1 = yt.plot_2d(ds, "Density",origin='native')
p1.annotate_grids()
p1.set_log("Density", False)
p1.show()

p2 = yt.plot_2d(ds, "Velocity_0",origin='native')
p2.annotate_grids()
p2.set_log("Velocity_0", False)
p2.show()

p3 = yt.plot_2d(ds, "Velocity_1",origin='native')
p3.annotate_grids()
p3.set_log("Velocity_1", False)
p3.show()

#p4 = yt.plot_2d(ds, "yvelnew",origin='native')
#p4.annotate_grids()
#p4.set_log("rhs", False)
#p4.show()

extent = [g.LeftEdge[0]-g.dds[0]/2, g.RightEdge[0]+g.dds[0]/2, g.LeftEdge[1]-g.dds[1]/2, g.RightEdge[1]+g.dds[1]/2]
plt.figure()
plt.clf()
plt.gca().set_aspect('equal')
plt.imshow(a_nd, origin='lower', extent=extent)
plt.colorbar()

