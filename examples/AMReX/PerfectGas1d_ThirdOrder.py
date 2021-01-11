import math

k0 = 0.017746/math.sqrt(300)*5000 #constant thermal conductivity
v0 = 2.23e-5/math.sqrt(293) *5000 #constant viscosity
kB = 1.38064852e-23 #Boltzmann constant
m = 6.6335209e-26 #gas(He) particles' mass
R = kB/m #gas constant
gamma = 1.6696 #specific heat ratio

Equation = {
  'gamma': gamma,
  'R_specific': R,
  'k_0': k0,
  'eta_0': v0,
}

GridGeometry = {
  'cell_dimensions': [200, 1, 1],
  'coordinates': {
    'lower': [0.0, 0.0, 0.0],
    'upper': [2.0, 2.0, 2.0]
  },
  'periodicity': [1, 1, 1]
}

PatchHierarchy = {
  'max_number_of_levels': 1,
  'refine_ratio': [2, 2, 2],
  'blocking_factor': [32, 32, 32],
  'max_grid_size': [256, 256, 256]
}

RunOptions = {
  'cfl': 0.5,
  'final_time': 1.0,
  'max_cycles': 2000
}

Output = {
  'outputs': [
    {
      'type': 'CounterOutput',
      'frequencies': [100]
    },
    {
      'type': 'HDF5',
      'path': 'WithDiffusion.h5',
      'frequencies': [1],
      'intervals': [1e-5]
    },
  ]
}
