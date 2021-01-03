import math

plenum_x_n_cells = 512
tube_blocking_factor = 8
plenum_blocking_factor = 8

n_level = 1

n_tubes = 6
r_tube = 0.015

D = 2.0 * r_tube

r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.2
r_tube_center = 2.0 * r_inner
alpha = 2.0 * math.pi / n_tubes

inlet_length = 3.0 * D # [m]

plenum_x_upper = 10.0 * D
plenum_x_lower = -inlet_length
plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_length = plenum_x_length # [m]
tube_length = 1.0 # [m]

plenum_max_grid_size = max(plenum_blocking_factor, 1024)

plenum_domain_length = plenum_length + inlet_length
tube_domain_length = tube_length - inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length

# plenum_yz_upper = +(r_outer + 0.01)
# plenum_yz_lower = -plenum_yz_upper

# plenum_yz_length = plenum_yz_upper - plenum_yz_lower


plenum_y_lower = - 0.0 * D
plenum_y_upper = + 15.0 * D
plenum_y_length = plenum_y_upper - plenum_y_lower

plenum_z_lower = plenum_y_lower
plenum_z_upper = plenum_y_upper
plenum_z_length = plenum_z_upper - plenum_z_lower

plenum_y_over_x_ratio = plenum_y_length / plenum_x_length

plenum_y_n_cells = plenum_x_n_cells * plenum_y_over_x_ratio
plenum_y_n_cells -= plenum_y_n_cells % plenum_blocking_factor
plenum_y_n_cells = int(plenum_y_n_cells)

plenum_z_over_x_ratio = plenum_z_length / plenum_x_length

plenum_z_n_cells = plenum_x_n_cells * plenum_z_over_x_ratio
plenum_z_n_cells -= plenum_z_n_cells % plenum_blocking_factor
plenum_z_n_cells = int(plenum_z_n_cells)

tube_n_cells = plenum_x_n_cells * tube_over_plenum_length_ratio
tube_n_cells -= tube_n_cells % tube_blocking_factor
tube_n_cells = int(tube_n_cells)

RunOptions = {
  'cfl': 0.8,
  'final_time': 10.0,
  'max_cycles': -1
}

# checkpoint = '/Users/maikel/Development/FiniteVolumeSolver/build_3d/MultiTube/Checkpoint/000000063'
checkpoint = ''

Plenum = {
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [plenum_x_n_cells, plenum_y_n_cells, plenum_z_n_cells],
    'coordinates': {
      'lower': [plenum_x_lower, plenum_y_lower, plenum_z_lower],
      'upper': [plenum_x_upper, plenum_y_upper, plenum_z_upper],
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [plenum_blocking_factor, plenum_blocking_factor, plenum_blocking_factor],
    'max_grid_size': [plenum_max_grid_size, plenum_max_grid_size, plenum_max_grid_size],
    'ngrow_eb_level_set': 9,
    'remove_covered_grids': False,
    'n_proper': 1,
    'n_error_buf': [0, 0, 0]
  },
  'InletGeometry': {
    'r_start': r_tube,
    'r_end': 2.0*r_tube
  },
  'BladeGeometry': {
    'dBox': 2.0 * r_tube,
    'dy': 2.0 * r_tube,
    'y0': 1000.0 * r_tube,
    'x0': 1000.0 * r_tube
  },
  'InitialCondition': {
    'left': {
      'density': 1.0 / 2.5,
      'temperature': 2.5
    },
    'right': {
      'pressure': 1.0 / 2.5,
      'temperature': 2.5
    },
  },
  'PressureBoundary': {
    'outer_pressure': 1.0
  }
}

def TubeCenterPoint(x0, k, alpha):
  return [x0, 0.0, 0.0]

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

Tube = {
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [tube_n_cells, 1, 1],
    'coordinates': {
      'lower': LowerX(-tube_length, 0, alpha),
      'upper': UpperX(-inlet_length, 0, alpha),
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [tube_blocking_factor, 1, 1],
    'refine_ratio': [2, 1, 1],
    'n_proper': 1,
    'n_error_buf': [4, 0, 0]
  }
}

def OuterProbe(x0, k, alpha):
  return [x0, r_outer - 0.002, r_tube_center * math.sin(k * alpha)]

Output = { 
  'outputs': [{
    'type': 'Plotfiles',
    'directory': 'ConvergentNozzle/Plotfiles/',
    'intervals': [0.01],
    #'frequencies': [1]
  }, {
    'type': 'Checkpoint',
    'directory': 'ConvergentNozzle/Checkpoint/',
    #'intervals': [1e-3],
    'frequencies': [100]
  }, {
   'type': 'CounterOutput',
   'intervals': [RunOptions['final_time'] / 3.0]
  }
  ]
}
