import math

RunOptions = {
  'cfl': 0.5,
  'final_time': 0.1 / math.sqrt(101325.),
  # 'final_time': 0.2,
}

GridGeometry = {
  'cell_dimensions': [128, 128, 1],
  'coordinates': {
    'lower': [-0.50, -0.50, -0.50],
    'upper': [+0.50, +0.50, +0.50],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
 'max_number_of_levels': 4,
 'refine_ratio': [2, 2, 1],
 'blocking_factor': [32, 32, 1],
 'max_grid_size': [32, 32, 1],
}

reconstruction = "Characteristics"
# reconstruction = "Primitive"
# reconstruction = "Conservative"
# reconstruction = "PerfectGas"

paths = {
  'Conservative': './SodShockTube/Conservative',
  'Primitive': './SodShockTube/Primitive',
  'Characteristics': './SodShockTube/Characteristics',
  'HLLE': './SodShockTube/HLLE'
}

Output = {
  'outputs': [{
    'type': 'Plotfile',
    'path': '{}.h5'.format(paths[reconstruction]),
    'directory': '{}_Plotfiles2'.format(paths[reconstruction]),
    'intervals': [RunOptions['final_time'] / 10.0],
    # 'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [5]
  }]
}
