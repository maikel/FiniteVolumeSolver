n_level = 1

blocking_factor = 8
max_grid_size = max(blocking_factor, 64)

Mach_number = 1.9

RunOptions = {
  'cfl': 0.4,
  'final_time': 0.0002,
  'max_cycles': -1,
  'do_backup': 0
}

n_cells = 128

GridGeometry = {
  'cell_dimensions': [n_cells, int(n_cells / 2), int(n_cells / 2)],
  'coordinates': {
    'lower': [-0.05, +0.000, 0.000],
    'upper': [+0.10, +0.075, 0.075],
  },
  'periodicity': [0, 0, 0]
}


InitialCondition = {
  'left': {
    'moles': 'O2:0.2, N2:0.8',
    'temperature': 2000,
    'pressure': 4.0 * 101325.0
  }, 'right': {
    'moles': 'O2:0.2, N2:0.8',
    'temperature': 300,
    'pressure': 101325.0
  }
}

PatchHierarchy = {
  'max_number_of_levels': n_level, 
  'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
  'max_grid_size': [max_grid_size, max_grid_size, max_grid_size],
  'ngrow_eb_level_set': 5,
  'remove_covered_grids': True
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'LinearShock3d/',
    'intervals': [1e-4],
    # 'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [25]
  }]
}
