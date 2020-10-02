import math

RunOptions = {
  'cfl': 0.5,
  'final_time': 0.001,
}

GridGeometry = {
  'cell_dimensions': [200, 1, 1],
  'coordinates': {
    'lower': [-1.50, -0.015, -0.015],
    'upper': [-0.00, +0.015, +0.015],
  },
  'periodicity': [1, 1, 1]
}

PatchHierarchy = {
 'max_number_of_levels': 1,
 'refine_ratio': [2, 1, 1],
 'blocking_factor': [8, 1, 1],
 'max_grid_size': [100, 1, 1],
}

reconstruction = "Primitive"

Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': './IdealGasMix_Prim.h5',
    'intervals': [5e-5],
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  }]
}
