RunOptions = {
  'cfl': 0.5,
  'final_time': 0.001,
  'max_cycles': 0
}

GridGeometry = {
  'cell_dimensions': [64, 64, 64],
  'coordinates': {
    'lower': [-0.06, -0.15, -0.15],
    'upper': [0.24, 0.15, 0.15],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': 2, 
  'blocking_factor': [8, 8, 8],
  'max_grid_size': [32, 32, 32],
  'ngrow_eb_level_set': 3
}

InitialCondition = {
  'left': {
    'pressure': 3.0 * 101325.0
  }, 'right': {
    'pressure': 1.0 * 101325.0
  }
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'LinearShock3d/',
    # 'intervals': [1e-4],
    'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [1]
  }]
}
