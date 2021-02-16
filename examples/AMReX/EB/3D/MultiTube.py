import math

plenum_x_n_cells = 64
tube_blocking_factor = 32
plenum_blocking_factor_x = 64
plenum_blocking_factor_y = 64
plenum_blocking_factor_z = 64

n_level = 1

n_tubes = 6
r_tube = 0.015
r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.389
r_tube_center = 2.0 * r_inner
alpha = 2.0 * math.pi / n_tubes

tube_length = 2.083 - 0.50 # [m]
inlet_length = 0.1 # [m]
plenum_length = 0.50 # [m]

plenum_max_grid_size_x = max(plenum_blocking_factor_x, 64)
plenum_max_grid_size_y = max(plenum_blocking_factor_y, 64)
plenum_max_grid_size_z = max(plenum_blocking_factor_z, 64)

plenum_domain_length = plenum_length + inlet_length
tube_domain_length = tube_length - inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length

# plenum_yz_upper = +(r_outer + 0.01)
# plenum_yz_lower = -plenum_yz_upper

# plenum_yz_length = plenum_yz_upper - plenum_yz_lower

plenum_x_upper = plenum_length + 0.1
plenum_x_lower = -inlet_length
plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_y_lower = - (r_outer + 0.02)
plenum_y_upper = + (r_outer + 0.02)
plenum_y_length = plenum_y_upper - plenum_y_lower

plenum_z_lower = - (r_outer + 0.02)
plenum_z_upper = + (r_outer + 0.02)
plenum_z_length = plenum_z_upper - plenum_z_lower

plenum_y_over_x_ratio = plenum_y_length / plenum_x_length

plenum_y_n_cells = plenum_x_n_cells * plenum_y_over_x_ratio
plenum_y_n_cells -= plenum_y_n_cells % plenum_blocking_factor_y
plenum_y_n_cells = int(plenum_y_n_cells)

plenum_z_over_x_ratio = plenum_z_length / plenum_x_length

plenum_z_n_cells = plenum_x_n_cells * plenum_z_over_x_ratio
plenum_z_n_cells -= plenum_z_n_cells % plenum_blocking_factor_z
plenum_z_n_cells = int(plenum_z_n_cells)

tube_n_cells = plenum_x_n_cells * tube_over_plenum_length_ratio
tube_n_cells -= tube_n_cells % tube_blocking_factor
tube_n_cells = int(tube_n_cells)

RunOptions = {
  'cfl': 0.8,
  'final_time': 0.009,
  'max_cycles': -1,
  'do_backup': 0,
}

# checkpoint = '/Users/maikel/Development/FiniteVolumeSolver/build_3d/MultiTube/Checkpoint/000000063'
# checkpoint = ''
checkpoint = '/srv/public/Maikel/FiniteVolumeSolver/build_3D-Release/MultiTube/Checkpoint/000000466'

Plenum = {
  'checkpoint': checkpoint,
  'InletGeometry': {
    'r_start': r_tube,
    'r_end': 0.0225
  },
  'BlockGeometry': {
    'factor': 0.0,
    'width': 10e-3
  },
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
    'blocking_factor': [plenum_blocking_factor_x, plenum_blocking_factor_y, plenum_blocking_factor_z],
    'max_grid_size': [plenum_max_grid_size_x, plenum_max_grid_size_y, plenum_max_grid_size_z],
    'ngrow_eb_level_set': 9,
    'cutcell_load_balance_weight': 2.,
    'remove_covered_grids': True,
  },
  'IntegratorContext': {
    'scratch_gcw': 2,
    'flux_gcw': 0,
    'regrid_frequency': 1,
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
  return [x0, r_tube_center * math.cos(k * alpha), r_tube_center * math.sin(k * alpha)]

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

fill_fraction = 1.0
measurement_position = - (1.0 - fill_fraction) * tube_length

ignition_position = -tube_length + 0.4
fuel_offsets = [0.005, 0.020, 0.005, 0.020, 0.005, 0.020]
ignition_offsets = [fuel_offset + 0.015 for fuel_offset in fuel_offsets]

IgniteDetonation = [{
  'interval': 0.030,
  'measurement_position': measurement_position,
  'equivalence_ratio_criterium': 0.9,
  'position': ignition_position,
  'offset': offset,
} for offset in ignition_offsets]

Tubes = [{
  'checkpoint': checkpoint,
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
  'IntegratorContext': {
    'scratch_gcw': 2,
    'flux_gcw': 0,
    'regrid_frequency': 1,
  },
  'PressureValveBoundary': {
    'prefix': 'PressureValve-{}'.format(i),
    'efficiency': 1.0,
    'open_at_interval': 0.030,
    'offset': fuel_offsets[i],
    'fuel_measurement_position': measurement_position,
    'fuel_measurement_criterium': 0.9,
    'pressure_value_which_opens_boundary': 101325.0,
    'pressure_value_which_closes_boundary': 3.0e5,
    'oxygen_measurement_position': -0.5,
    'oxygen_measurement_criterium': 0.1,
    'equivalence_ratio': 1.0,
    'massflow_boundary': {
      'coarse_inner_box': { 
        'lower': [0, 0, 0], 
        'upper': [1, 0, 0] 
      },
      'side': 0,
      'direction': 0,
      'required_massflow': 900.0 / 6.0 / 3600.0,
      'surface_area': math.pi * r_tube * r_tube
    }
  }
} for i in range(0, n_tubes)]

def OuterProbe(x0, k, alpha):
  return [x0, (r_outer - 0.002) * math.cos(k *alpha), r_tube_center * math.sin(k * alpha)]

slices = [{
    'type': 'HDF5',
    'which_block': i + 1,
    'path': 'MultiTube/Slices/Tube_{}.h5'.format(i),
    'intervals': [1e-4],
    'frequencies': range(450, 500)
  } for i in range(0, 6)]

Output = { 
  'outputs': [{ 
    'type': 'Plotfile',
    'directory': 'MultiTube/Plenum/',
    'intervals': [5e-4],
    # 'frequencies': [1]
  },
  # , {
    # 'type': 'HDF5',
    # 'which_block': 0,
    # 'path': 'MultiTube/Slices/Plenum.h5',
    # 'intervals': [1e-4],
    # 'frequencies': [1],
    # 'box': {
      # 'lower': [0, 0, int(plenum_z_n_cells / 2)],
      # 'upper': [plenum_x_n_cells - 1, plenum_y_n_cells - 1, int(plenum_z_n_cells / 2)]
    # }
  # },
  {
    'type': 'LogProbes',
    'directory': 'MultiTube/Probes/',
    'frequencies': [10],
    'Plenum': {
      'filename': 'MultiTube/Probes/Plenum.h5',
      'coordinates': [
        OuterProbe(0.07, 0, alpha), OuterProbe(0.16, 0, alpha), OuterProbe(0.25, 0, alpha), OuterProbe(0.34, 0, alpha), OuterProbe(0.43, 0, alpha),
        OuterProbe(0.07, 1, alpha), OuterProbe(0.16, 1, alpha), OuterProbe(0.25, 1, alpha), OuterProbe(0.34, 1, alpha), OuterProbe(0.43, 1, alpha),
        OuterProbe(0.07, 2, alpha), OuterProbe(0.16, 2, alpha), OuterProbe(0.25, 2, alpha), OuterProbe(0.34, 2, alpha), OuterProbe(0.43, 2, alpha),
        OuterProbe(0.07, 3, alpha), OuterProbe(0.16, 3, alpha), OuterProbe(0.25, 3, alpha), OuterProbe(0.34, 3, alpha), OuterProbe(0.43, 3, alpha),
        OuterProbe(0.07, 4, alpha), OuterProbe(0.16, 4, alpha), OuterProbe(0.25, 4, alpha), OuterProbe(0.34, 4, alpha), OuterProbe(0.43, 4, alpha),
        OuterProbe(0.07, 5, alpha), OuterProbe(0.16, 5, alpha), OuterProbe(0.25, 5, alpha), OuterProbe(0.34, 5, alpha), OuterProbe(0.43, 5, alpha),
        TubeCenterPoint(-0.07, 0, alpha),  TubeCenterPoint(-0.07, 1, alpha), TubeCenterPoint(-0.07, 2, alpha), TubeCenterPoint(-0.07, 3, alpha), TubeCenterPoint(-0.07, 4, alpha), TubeCenterPoint(-0.07, 5, alpha)
      ],
    },
    'Tube': {
      'n_tubes': 6,
      'filename': 'MultiTube/Probes/Tube',
      'coordinates': [
        TubeCenterPoint(-0.9, 0, alpha), TubeCenterPoint(-0.45, 0, alpha), TubeCenterPoint(-0.3, 0, alpha), TubeCenterPoint(-0.15, 0, alpha), 
        TubeCenterPoint(-0.9, 1, alpha), TubeCenterPoint(-0.45, 1, alpha), TubeCenterPoint(-0.3, 1, alpha), TubeCenterPoint(-0.15, 1, alpha),
        TubeCenterPoint(-0.9, 2, alpha), TubeCenterPoint(-0.45, 2, alpha), TubeCenterPoint(-0.3, 2, alpha), TubeCenterPoint(-0.15, 2, alpha),
        TubeCenterPoint(-0.9, 3, alpha), TubeCenterPoint(-0.45, 3, alpha), TubeCenterPoint(-0.3, 3, alpha), TubeCenterPoint(-0.15, 3, alpha),
        TubeCenterPoint(-0.9, 4, alpha), TubeCenterPoint(-0.45, 4, alpha), TubeCenterPoint(-0.3, 4, alpha), TubeCenterPoint(-0.15, 4, alpha),
        TubeCenterPoint(-0.9, 5, alpha), TubeCenterPoint(-0.45, 5, alpha), TubeCenterPoint(-0.3, 5, alpha), TubeCenterPoint(-0.15, 5, alpha),
      ],
    }
  }, {
    'type': 'Checkpoint',
    'directory': 'MultiTube/Checkpoint/',
    'intervals': [1e-3],
    'frequencies': range(450, 500)
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  }]
}

Output['outputs'].extend(slices)
