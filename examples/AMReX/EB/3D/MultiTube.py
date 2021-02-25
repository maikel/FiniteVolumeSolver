import math

print(args)
if 'args' in globals() and len(args) > 0:
  ignition_position = float(args[0])
else:
  ignition_position = -tube_length + 0.4

if 'args' in globals() and len(args) > 1:
  checkpoint = args[1]
else:
  checkpoint = ''

print(args)
print(ignition_position)

plenum_x_n_cells = 64
tube_blocking_factor = 32
plenum_blocking_factor_x = 16
plenum_blocking_factor_y = 16
plenum_blocking_factor_z = 16

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
  'cfl': 0.7,
  'final_time': 0.040,
  'max_cycles': -1,
  'do_backup': 1,
}

# checkpoint = '/Users/maikel/Development/FiniteVolumeSolver/build_3d/MultiTube/Checkpoint/000000063'
# checkpoint = ''
# checkpoint = '/srv/public/Maikel/FiniteVolumeSolver/build_3D-Release/MultiTube/Checkpoint/000001166'

Plenum = {
  'checkpoint': checkpoint,
  'InletGeometry': {
    'r_start': r_tube,
    'r_end': 1.5 * r_tube
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

ignition_position = -tube_length + 0.2
fuel_offsets = [0.004, 0.050, 0.05, 0.050, 0.05, 0.05]
ignition_offsets = [fuel_offset + 0.030 for fuel_offset in fuel_offsets]

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
      'required_massflow': 0.2 / 6.0,
      'surface_area': math.pi * r_tube * r_tube
    }
  }
} for i in range(0, n_tubes)]

def OuterProbe(x0, k, alpha):
  return [x0, (r_outer - 0.002) * math.cos(k *alpha), r_tube_center * math.sin(k * alpha)]

slices = [{
    'type': 'HDF5',
    'which_block': i + 1,
    'path': 'MultiTube2/Slices/Tube_{}.h5'.format(i),
    'intervals': [1e-4],
    #'frequencies': range(450, 500)
  } for i in range(0, 6)]

import itertools
xvalues = [0.07 + 0.9 * i for i in range(0, 6)]
probes = [OuterProbe(xvalue, tube, alpha) for (xvalue, tube) in itertools.product(xvalues, range(0, 6))]
probes.extend([TubeCenterPoint(-0.07, tube, alpha) for tube in range(0, 6)])

txvalues = [-0.15 - 0.15 * i for i in range(0,7)]
tprobes = [TubeCenterPoint(xvalue, tube, alpha) for (xvalue, tube) in itertools.product(txvalues, range(0, 6))]

Output = { 
  'outputs': [{ 
    'type': 'Plotfile',
    'directory': 'MultiTube2/Plenum/',
    'intervals': [],
    #'frequencies': range(450, 500)
  },
  # , {
    # 'type': 'HDF5',
    # 'which_block': 0,
    # 'path': 'MultiTube2/Slices/Plenum.h5',
    # 'intervals': [1e-4],
    # 'frequencies': [1],
    # 'box': {
      # 'lower': [0, 0, int(plenum_z_n_cells / 2)],
      # 'upper': [plenum_x_n_cells - 1, plenum_y_n_cells - 1, int(plenum_z_n_cells / 2)]
    # }
  # },
  {
    'type': 'LogProbes',
    'directory': 'MultiTube2/Probes/',
    'frequencies': [10],
    'Plenum': {
      'filename': 'MultiTube2/Probes/Plenum.h5',
      'coordinates': probes,
    },
    'Tube': {
      'n_tubes': 6,
      'filename': 'MultiTube2/Probes/Tube',
      'coordinates': tprobes,
    }
  }, {
    'type': 'Checkpoint',
    'directory': 'MultiTube2/Checkpoint/',
    'intervals': [1e-3],
    #'frequencies': range(450, 500)
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  }]
}

Output['outputs'].extend(slices)
