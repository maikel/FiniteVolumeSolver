import math

RunOptions = {
  'cfl': 0.8,
  # 'final_time': 0.2 / math.sqrt(101325.),
  'final_time': 5.0,
  # 'max_cycles': 1,
}

dx = 5e-3
x_len = 1.0
n_cells = int(x_len / dx)

GridGeometry = {
  'cell_dimensions': [n_cells, 1, 1],
  'coordinates': {
    'lower': [-0.0, -0.015, -0.015],
    'upper': [+1.00, +0.015, +0.015],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
 'max_number_of_levels': 1,
 'refine_ratio': [2, 1, 1],
 'blocking_factor': [8, 1, 1],
 'max_grid_size': [n_cells, 1, 1],
}

reconstruction = "Characteristics"
# reconstruction = "Primitive"
# reconstruction = "NoReconstruct"

paths = {
  'HLLE': './HLLE.h5',
  'NoReconstruct': './HLLEM.h5',
  'Conservative': './Conservative.h5',
  'ConservativeM': './ConservativeM.h5',
  'Primitive': './Primitive.h5',
  'Characteristics': './Characteristics.h5',
  'PerfectGas': './PerfectGas1d_cfl05.h5'
}

Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': paths[reconstruction],
    'intervals': [0.01],
    #'frequencies': [1],
  },
  ]
}
