import math

RunOptions = {
  'cfl': 0.9,
  'final_time': 0.2# / math.sqrt(1.),
  # 'final_time': 1.0,
  # 'max_cycles': 1,
}

dx = 0.5e-2
x_len = 1.0
n_cells = int(x_len / dx)

GridGeometry = {
  'cell_dimensions': [n_cells, 1, 1],
  'coordinates': {
    'lower': [-0.50, -0.015, -0.015],
    'upper': [+0.50, +0.015, +0.015],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
 'max_number_of_levels': 1,
 'refine_ratio': [2, 1, 1],
 'blocking_factor': [8, 1, 1],
 'max_grid_size': [n_cells, 1, 1],
}

FluxMethod = {
  'reconstruction': 'Characteristics',
  'limiter': 'VanLeer',
  'base_method': 'HLLE'
}

Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': 'PerfectGas1d_Char.h5',
    'intervals': [RunOptions['final_time'] / 10.],
    #'frequencies': [1],
  },{
    'type': 'CounterOutput',
    'intervals': [RunOptions['final_time']],
    #'frequencies': [1],
  },
  ]
}
