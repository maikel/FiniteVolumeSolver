import math

RunOptions = {
  'cfl': 0.8,
  'final_time': 0.001,
}

GridGeometry = {
  'cell_dimensions': [800, 1, 1],
  'coordinates': {
    'lower': [-1.00, -0.015, -0.015],
    'upper': [+1.00, +0.015, +0.015],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
 'max_number_of_levels': 1,
 'refine_ratio': [2, 1, 1],
 'blocking_factor': [8, 1, 1],
 'max_grid_size': [100, 1, 1],
}

reconstruction = "Characteristics"
# reconstruction = "Primitive"
# reconstruction = "Conservative"

Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': './IdealGasMix_Char.h5',
    # 'path': './IdealGasMix_Prim.h5',
    # 'path': './IdealGasMix_Cons.h5',
    'intervals': [5e-5],
  }]
}
