RunOptions = {
  'cfl': 0.5,
  'final_time': 0.02,
  'max_cycles': 10
}

GridGeometry = {
  'cell_dimensions': [32, 32, 32],
  'coordinates': {
    'lower': [-0.03, -0.15, -0.15],
    'upper': [0.27, 0.15, 0.15],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': 1, 
  'blocking_factor': [16, 16, 16],
  'max_grid_size': [32, 32, 32],
  'ngrow_eb_level_set': 5
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'MultiTube/',
    # 'intervals': [1e-4],
    'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [1]
  }]
}
