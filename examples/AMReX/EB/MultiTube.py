import math

n_tubes = 6
r_tube = 0.015
r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.385
r_tube_center = 0.5 * r_inner + 0.5 * r_outer
alpha = 2.0 * math.pi / n_tubes

RunOptions = {
  'cfl': 0.95,
  'final_time': 1.0,
  'max_cycles': 0
}

Plenum = {
  'GridGeometry': {
    'cell_dimensions': [64, 64, 64],
    'coordinates': {
      'lower': [-0.03, -0.53 / 2, -0.53 / 2],
      'upper': [ 0.50, +0.53 / 2, +0.53 / 2],
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': 1, 
    'blocking_factor': [64, 64, 64],
    'max_grid_size': [64, 64, 64]
  },
  'IsentropicPressureBoundary': {
    'outer_pressure': 101325.0,
    'coarse_inner_box': { 'lower': [31, 0, 0], 'upper': [31, 31, 31] },
    'side': 1,
    'direction': 0
  },
  'InitialCondition': {
    'moles': 'O2:20,N2:80',
    'temperature': 300.0,
    'pressure': 101325.0
  }
}

def TubeCenterPoint(x0, k, alpha):
  return [x0, r_tube_center * math.sin(k * alpha), r_tube_center * math.cos(k * alpha)]

def LowerX(x0, k, alpha):
  center = TubeCenterPoint(x0, k, alpha)
  center[1] -= r_tube
  center[2] -= r_tube
  return center

def UpperX(x0, k, alpha):
  center = TubeCenterPoint(x0, k, alpha)
  center[1] += r_tube
  center[2] += r_tube
  return center

Tubes = [{
  'GridGeometry': {
    'cell_dimensions': [64, 1, 1],
    'coordinates': {
      'lower': LowerX(-0.56, i, alpha),
      'upper': UpperX(-0.03, i, alpha),
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': 1, 
    'blocking_factor': [32, 1, 1]
  },
  'PressureValveBoundary': {
    'efficiency': 1.0,
    'open_at_interval': 0.03333333,
    'offset': 0.0,
    'fuel_measurement_position': -0.2,
    'fuel_measurement_criterium': 0.8,
    'pressure_value_which_opens_boundary': 101325.0,
    'pressure_value_which_closes_boundary': 3.0e5,
    'oxygen_measurement_position': -0.5,
    'oxygen_measurement_criterium': 0.1,
    'equivalence_ratio': 1.0
  }
} for i in range(0, n_tubes)]


Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'MultiTube/Test/',
    'intervals': [1e-4],
  }]
}

IgniteDetonation = {
  'interval': 0.03333333,
  'measurement_position': -0.2,
  'equivalence_ratio_criterium': 0.8,
  'position': 0.8
}


