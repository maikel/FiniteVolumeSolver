RunOptions = {
  'cfl': 0.4,# should be between in (0, 1]
  'final_time': 0.002,
  'max_cycles': -1, # -1 means infinite and 0 means only initial condition
}

n_cells = 256
n_levels = 1

GridGeometry = {
  'cell_dimensions': [n_cells, int(n_cells / 2), 1],
  'coordinates': {
    'lower': [-0.05, +0.000, -1.0],
    'upper': [+0.10, +0.075, +1.0],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_levels, 
  'n_error_buf': [0, 0, 0],
  'ngrow_eb_level_set': 9,
  'remove_covered_grids': False,
  'n_proper': 1,
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'ReferenceData/LinearShockAxi_{}x{}/'.format(n_cells, n_levels),
    'intervals': [0.0001],
  }, {
    'type': 'HDF5',
    'path': 'ReferenceData/LinearShockAxi_{}x{}.h5'.format(n_cells, n_levels),
    'intervals': [0.00001],
  }, {
    'type': 'CounterOutput',
    'intervals': [0.0001]
  }]
}
