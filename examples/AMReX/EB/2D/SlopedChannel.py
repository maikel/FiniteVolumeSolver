import math

RunOptions = {
  'cfl': 0.4,
  'final_time': 0.0015,
  'max_cycles': -1, # means infinite and 0 means only initial condition
}

def AlignForBlockingFactor(n, blocking_factor):
  return n - n % blocking_factor

blocking_factor = 1

factor=8
n_cells_x = 50*factor
n_cells_y = AlignForBlockingFactor(50*factor, blocking_factor)
n_levels = 1

GridGeometry = {
  'cell_dimensions': [n_cells_x, n_cells_y, 1],
  'coordinates': {
    'lower': [+0.00, +0.00, +0.00],
    'upper': [+0.10, +0.10, +0.07],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_levels, 
  'n_error_buf': [1, 1, 1],
  'max_grid_size': [1024, 1024, 1024],
  'ngrow_eb_level_set': 9,
  'remove_covered_grids': False,
  'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
  'n_proper': 1,
}

# reconstruction = "HLLE"
reconstruction = "Characteristics"
# reconstruction = "Conservative"

origin = 0.035
theta = math.pi * 45.0 / 180.0

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'ReferenceData/SlopedChannel_old_{}_{}x{}-{}/'.format(reconstruction, n_cells_x, n_cells_y, n_levels),
    'intervals': [1e-4],
    # 'frequencies': [1] 
  },
  {
    'type': 'HDF5',
    'path': 'ReferenceData/SlopedChannel_old_{}_{}x{}-{}.h5'.format(reconstruction, n_cells_x, n_cells_y, n_levels),
    'intervals': [1e-4],
    # 'frequencies': [1] 
  }]
}
