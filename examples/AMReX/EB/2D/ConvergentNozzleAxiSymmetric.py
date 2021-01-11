import math

plenum_x_n_cells = 256
tube_blocking_factor = 8
plenum_blocking_factor = 8

n_level = 2

n_tubes = 6
r_tube = 0.015

D = 2.0 * r_tube

r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.2
r_tube_center = 2.0 * r_inner
alpha = 2.0 * math.pi / n_tubes

inlet_length = 3.0 * D # [m]

plenum_x_upper = 20.0 * D
plenum_x_lower = -inlet_length
plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_length = plenum_x_length # [m]
tube_length = 2.083 - 0.5 # [m]

plenum_max_grid_size = max(plenum_blocking_factor, 64)

plenum_domain_length = plenum_length + inlet_length
tube_domain_length = tube_length - inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length

# plenum_yz_upper = +(r_outer + 0.01)
# plenum_yz_lower = -plenum_yz_upper

# plenum_yz_length = plenum_yz_upper - plenum_yz_lower


plenum_y_lower = 0.0
plenum_y_upper = + 10.0 * D
plenum_y_length = plenum_y_upper - plenum_y_lower

plenum_z_lower = - 10.0 * D
plenum_z_upper = + 10.0 * D
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
  'cfl': 0.4,
  'final_time': 0.021,
  'max_cycles': -1
}

#checkpoint = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Release/ConvergentNozzleAxi/Checkpoint/000012985'
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
    'ngrow_eb_level_set': 5,
    'remove_covered_grids': False,
    'n_proper': 1,
    'n_error_buf': [0, 0, 0]
  },
  'InletGeometry': {
    'r_start': r_tube,
    'r_end': r_tube
  },
  'IsentropicPressureBoundary': {
    'outer_pressure': 101325.0,
    'coarse_inner_box': { 
      'lower': [plenum_x_n_cells - 1, 0, 0],
      'upper': [plenum_x_n_cells - 1, plenum_y_n_cells - 1, plenum_z_n_cells - 1] 
    },
    'side': 1,
    'direction': 0
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
  },
  'PressureValveBoundary': {
    'prefix': 'PressureValve',
    'efficiency': 1.0,
    'open_at_interval': 0.03333333,
    'offset': 0.01,
    'fuel_measurement_position': -0.15,
    'fuel_measurement_criterium': 0.9,
    'pressure_value_which_opens_boundary': 101325.0,
    'pressure_value_which_closes_boundary': 3.0e5,
    'equivalence_ratio': 1.0,
    'massflow_boundary': {
      'coarse_inner_box': { 
        'lower': [0, 0, 0], 
        'upper': [1, 0, 0] 
      },
      'side': 0,
      'direction': 0,
      'required_massflow': 100.0 / 3600.0, # kg / s
      'surface_area': math.pi * r_tube * r_tube
    }
  }
}

def OuterProbe(x0, k, alpha):
  return [x0, r_outer - 0.002, r_tube_center * math.sin(k * alpha)]

Output = { 
  'outputs': [
  {
    'type': 'HDF5',
    'which_block': 0,
    'path': 'ConvergentNozzleAxi/Plenum.h5',
    'intervals': [1e-4]
  },
  {
    'type': 'Plotfiles',
    'directory': 'ConvergentNozzleAxi/Plotfiles/',
    'intervals': [1e-4],
    # 'frequencies': [1]
  }, {
    'type': 'Checkpoint',
    'directory': 'ConvergentNozzleAxi2/Checkpoint/',
    'intervals': [1e-3],
    'frequencies': []
  }, {
    'type': 'CounterOutput',
    'frequencies': [1000]
  }]
}

IgniteDetonation = {
  'interval': 0.06,
  'offset': 0.018,
  'measurement_position': -0.3, # -0.45
  'equivalence_ratio_criterium': 0.9,
  'position': -0.8
}
