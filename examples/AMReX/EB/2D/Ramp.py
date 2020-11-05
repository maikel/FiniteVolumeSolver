import math

isRef = 0
factor = 1

final_time = 0.5

max_cycles = [1, 0]
mode = ['new', 'ref']

RunOptions = {
  'cfl': 0.4,
  'final_time': final_time,
  'max_cycles': 1 # -1 means infinite and 0 means only initial condition
}

def AlignForBlockingFactor(n, blocking_factor):
  return n - n % blocking_factor

blocking_factor = 1

n_cells_x = 60*factor
n_cells_y = AlignForBlockingFactor(40*factor, blocking_factor)
n_levels = 1

GridGeometry = {
  'cell_dimensions': [n_cells_x, n_cells_y, 1],
  'coordinates': {
    'lower': [-1.00, -1.00, -1.00],
    'upper': [+2.00, +1.00, +1.00],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_levels, 
  'n_error_buf': [1, 1, 1],
  'max_grid_size': [n_cells_x, n_cells_y, n_cells_x],
  'ngrow_eb_level_set': 9,
  'remove_covered_grids': False,
  'blocking_factor': [n_cells_x, n_cells_y, 1],
  'n_proper': 1,
}

Output = { 
  'outputs': [
    {
    'type': 'Plotfile',
    'directory': 'Debug_Ramp_{}x{}-{}/Plot/'.format(n_cells_x, n_cells_y, n_levels),
    'intervals': [final_time],
    # 'frequencies': [1] 
  },
  {
    'type': 'DebugOutput',
    'directory': 'Debug_Ramp_{}x{}-{}/'.format(n_cells_x, n_cells_y, n_levels),
    'intervals': [1e-4],
    'frequencies': [1] 
  },
  {
    'type': 'HDF5',
    'path': 'ReferenceData/Ramp_{}x{}-{}.h5'.format(n_cells_x, n_cells_y, n_levels),
    'intervals': [final_time],
    # 'frequencies': [1] 
  }]
}
