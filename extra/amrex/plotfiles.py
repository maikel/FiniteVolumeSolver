import yt
import numpy as np

def load(path_to_plotfile, vars, mask_boundary_cells=True, buffer_size=None, vfrac_cutoff=1e-15):
  ds = yt.load(path_to_plotfile)
  plot = yt.plot_2d(ds, vars, origin='native')
  if buffer_size:
    plot.set_buff_size(buffer_size); 
  else:
    plot.set_buff_size(ds.domain_dimensions)
  current_time = float(ds.current_time)
  bounds = np.array(plot.bounds)
  datas = [np.squeeze(np.array(plot.frb[var])) for var in vars]
  if mask_boundary_cells:
    vfrac = np.array(plot.frb["vfrac"])
    vfrac = np.squeeze(vfrac)
    datas = [np.ma.masked_where( vfrac<=vfrac_cutoff, data) for data in datas]
  return datas, current_time, bounds