import math

plenum_x_n_cells = 16 * 8
tube_blocking_factor = 8
plenum_blocking_factor = 8

n_level = 1

n_tubes = 6
r_tube = 0.015
r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.385
r_tube_center = 0.5 * r_inner + 0.5 * r_outer
alpha = 2.0 * math.pi / n_tubes

tube_length = 1.50 # [m]
inlet_length = 0.03 # [m]
plenum_length = 0.50 # [m]

plenum_max_grid_size = max(plenum_blocking_factor, 32)

plenum_domain_length = plenum_length + inlet_length
tube_domain_length = tube_length - inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length

# plenum_yz_upper = +(r_outer + 0.01)
# plenum_yz_lower = -plenum_yz_upper

# plenum_yz_length = plenum_yz_upper - plenum_yz_lower

plenum_x_upper = plenum_length
plenum_x_lower = -inlet_length

plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_yz_length = plenum_x_length
plenum_yz_lower = -plenum_yz_length * 0.5
plenum_yz_upper = +plenum_yz_length * 0.5

plenum_yz_over_x_ratio = plenum_yz_length / plenum_x_length

plenum_yz_n_cells = plenum_x_n_cells * plenum_yz_over_x_ratio
plenum_yz_n_cells -= plenum_yz_n_cells % plenum_blocking_factor
plenum_yz_n_cells = int(plenum_yz_n_cells)

tube_n_cells = plenum_x_n_cells * tube_over_plenum_length_ratio
tube_n_cells -= tube_n_cells % tube_blocking_factor
tube_n_cells = int(tube_n_cells)

RunOptions = {
  'cfl': 0.5,
  'final_time': 0.02,
  # 'max_cycles': 0
}

Plenum = {
  'GridGeometry': {
    'cell_dimensions': [plenum_x_n_cells, plenum_yz_n_cells, plenum_yz_n_cells],
    'coordinates': {
      'lower': [plenum_x_lower, plenum_yz_lower, plenum_yz_lower],
      'upper': [plenum_x_upper, plenum_yz_upper, plenum_yz_upper],
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [plenum_blocking_factor, plenum_blocking_factor, plenum_blocking_factor],
    'max_grid_size': [plenum_max_grid_size, plenum_max_grid_size, plenum_max_grid_size]
  },
  'IsentropicPressureBoundary': {
    'outer_pressure': 101325.0,
    'coarse_inner_box': { 
      'lower': [plenum_x_n_cells - 1, 0, 0],
      'upper': [plenum_x_n_cells - 1, plenum_yz_n_cells - 1, plenum_yz_n_cells - 1] 
    },
    'side': 1,
    'direction': 0
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
    'cell_dimensions': [tube_n_cells, 1, 1],
    'coordinates': {
      'lower': LowerX(-tube_length, i, alpha),
      'upper': UpperX(-inlet_length, i, alpha),
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [tube_blocking_factor, 1, 1],
    'refine_ratio': [2, 1, 1]
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
    # 'intervals': [1e-4],
    'frequencies': [1],
  }, {
    'type': 'LogProbes',
    'directory': 'MultiTube/Probes/',
    'frequencies': [1],
    'Plenum': {
      'filename': 'MultiTube/Probes/Plenum.h5',
      'coordinates': [TubeCenterPoint(0.0, 0, alpha)],
    },
    'Tube': {
      'filename': 'MultiTube/Probes/Tubes.h5',
      'coordinates': [TubeCenterPoint(-1.0, 0, alpha), TubeCenterPoint(-0.2, 0, alpha)],
    }
  }]
}

IgniteDetonation = {
  'interval': 0.03333333,
  'measurement_position': -0.2,
  'equivalence_ratio_criterium': 0.8,
  'position': -0.8
}


