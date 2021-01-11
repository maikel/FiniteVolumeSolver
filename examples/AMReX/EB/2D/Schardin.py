RunOptions = {
  'cfl': 0.4,
  'final_time': 0.0002,
  'max_cycles': -1, # -1 means infinite and 0 means only initial condition
}

def AlignForBlockingFactor(n, blocking_factor):
  return n - n % blocking_factor

blocking_factor = 4

n_cells_x = 200
n_cells_y = AlignForBlockingFactor(140, blocking_factor)
n_levels = 2

GridGeometry = {
  'cell_dimensions': [n_cells_x, n_cells_y, 1],
  'coordinates': {
    'lower': [+0.00, +0.00, +0.00],
    'upper': [+0.15, +0.10, +0.10],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_levels, 
  'n_error_buf': [0, 0, 0],
  'max_grid_size': [1024, 1024, 1024],
  'ngrow_eb_level_set': 9,
  'remove_covered_grids': False,
  'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
  'n_proper': 1,
}

reconstruction = "HLLEM"

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'ReferenceData/Schardin_{}x{}-{}/'.format(n_cells_x, n_cells_y, n_levels),
    'intervals': [0.00005],
  }, {
    'type': 'HDF5',
    'path': 'ReferenceData/Schardin_{}x{}-{}.h5'.format(n_cells_x, n_cells_y, n_levels),
    'intervals': [0.000005],
  }]
}
