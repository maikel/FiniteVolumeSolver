import numpy as np
import yt

"""
select yt log level:


example:
  import yt
  yt.funcs.mylog.setLevel(level)

  --> level: int or str
        Possible values by increasing level:
        0 or "notset"
        1 or "all"
        10 or "debug"
        20 or "info"
        30 or "warning"
        40 or "error"
        50 or "critical"
"""

"""
Creating derived data fields see: 
https://yt-project.org/doc/developing/creating_derived_fields.html#creating-derived-fields

example for 2D-velocityMagnitude:

  import yt
  def _velMag(field, data):
    velX = data["Momentum_0"]/data["Density"]
    velY = data["Momentum_1"]/data["Density"]

    magVel = np.sqrt( np.power(velX, 2) + np.power(velY, 2) )
    return magVel

  yt.add_field(
    name=('velocityMagnitude'),
    function=_velMag,
    take_log=False
  )

"""

def yt_load(path_to_plotfile, vars, mask_boundary_cells=True, buffer_size=None, vfrac_cutoff=1e-15):
  """
  Load the with 'vars' specified 2D data from the plt-files.

  Parameters
  ----------------------------------------
    path_to_plotfile:     string
                          path to the plt-files
    vars:                 list of strings
                          list of variables to extract data
    mask_boundary_cells:  boolean, optional
                          mask the cutcell data with the volume fraction
    buffer_size:          optional
                          buffers loading data
    vfrac_cutoff:         float, optional
                          the cutoff parameter for masking the cutcell data

  Returns
  ----------------------------------------
    datas:        numpy array
                  the loaded datas from the plt-file
    current_time: float
                  the time point from the data
    bounds:       list of floats
                  the physical extension from the dataset
  """
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

def yt_collect_child_grids(startGrid):
  """
  Extract all child grids from given yt base grid.
  See also: https://yt-project.org/doc/examining/low_level_inspection.html

  Parameters
  ----------------------------------------
    startGrid:            index.grid object from yt dataset object
                          the grid where to extract all child grids

  Returns
  ----------------------------------------
    grids:        python list
                  all extracted child grids
  """
  stack = [startGrid]
  grids = []
  while stack:
      currentGrid = stack.pop()
      grids.append(currentGrid)
      stack.extend(grid for grid in currentGrid.Children)
  return grids