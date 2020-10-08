import math

RunOptions = {
  'cfl': 0.5,
  'final_time': 0.1 / math.sqrt(101325.),
  # 'final_time': 0.2,
}

GridGeometry = {
  'cell_dimensions': [1024, 1024, 1],
  'coordinates': {
    'lower': [-0.50, -0.50, -0.50],
    'upper': [+0.50, +0.50, +0.50],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
 'max_number_of_levels': 1,
 'refine_ratio': [2, 2, 1],
 'blocking_factor': [32, 32, 1],
 'max_grid_size': [32, 32, 1],
}

# reconstruction = "HLLE"
# reconstruction = "Primitive"
# reconstruction = "Conservative"
# reconstruction = "ConservativeM"
reconstruction = "Characteristics"
# reconstruction = "PerfectGas"

paths = {
  'Conservative': './SodShockTube/Conservative',
  'ConservativeM': './SodShockTube/ConservativeM',
  'Primitive': './SodShockTube/Primitive',
  'Characteristics': './SodShockTube/Characteristics',
  'HLLE': './SodShockTube/HLLE',
  'HLLEM': './SodShockTube/HLLEM'
}

Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': '{}.h5'.format(paths[reconstruction]),
    'directory': '{}_Plotfiles2'.format(paths[reconstruction]),
    'intervals': [RunOptions['final_time'] / 10.0],
    # 'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [5]
  }]
}
