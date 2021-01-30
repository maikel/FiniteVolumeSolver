import math

#dx = 256
x_len = 1.0
n_cells = 256# int(x_len / dx)

RunOptions = {
  'cfl': 0.5 * 0.9 / float(n_cells / 64),
  # 'final_time': 0.2 / math.sqrt(101325.),
  'final_time': 5.0,
  # 'max_cycles': 1,
}

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

# ArrheniusKinetics = {
#   'Q': 8.0,
#   'EA': 11.0,
#   'B': 0.05,
#   'T_switch':1.10,
# }

CompressorState = {
  'pressure': 2.0,
  'temperature': 1.0
}

def Area(xi):
  A0  = 1.0
  A1  = 4.0 # Reference: 3.0; best: 4.0 
  xi0 = 0.0625 - 1.0
  xi1 = 0.75 - 1.0   # Reference: 0.5; best: 0.75
  Ax = 1.0 if xi < xi0 else A0 + (A1-A0)*(xi-xi0)/(xi1-xi0) if xi < xi1 else A1
  return Ax

FluxMethod = {
  'reconstruction': 'Characteristics',
  'axial_variation': Area
}


Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': 'Deflagration.h5',
    'intervals': [0.005],
    #'frequencies': [1],
  },{
    'type': 'CounterOutput',
    'intervals': [1.0],
    #'frequencies': [1],
  },
  ]
}
