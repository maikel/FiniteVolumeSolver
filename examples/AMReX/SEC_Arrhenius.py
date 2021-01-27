import math

RunOptions = {
  'cfl': 0.8,
  # 'final_time': 0.2 / math.sqrt(101325.),
  'final_time': 50.0,
  # 'max_cycles': 1,
}

dx = 1e-3
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

ArrheniusKinetics = {
  'Q': 8.0,
  'EA': 11.0,
  'B': 0.05,
  'T_switch':1.10,
}

reconstruction = "Characteristics"
#reconstruction = "Conservative"
#reconstruction = "Primitive"
#reconstruction = "NoReconstruct"

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
  },{
    'type': 'CounterOutput',
    'intervals': [1.0],
    #'frequencies': [1],
  },
  ]
}
