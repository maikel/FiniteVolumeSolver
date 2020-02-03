RunOptions = {
  'cfl': 0.5,
  'final_time': 0.001,
  'max_cycles': 0
}

CartesianGridGeometry = {
  'cell_dimensions': [200, 1],
  'coordinates': {
    'lower': [-1.50, -0.015],
    'upper': [-0.06, +0.015],
  },
  'periodicity': [0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': 1, 
  'refine_ratio': [2, 1],
  'blocking_factor': [8, 1],
  'max_grid_size': [120, 1],
  'periodicity': [0, 0]
}

TemperatureRamp = {
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
