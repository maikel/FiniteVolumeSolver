import math

RunOptions = {
  'cfl': 0.5,
  'final_time': 0.04,
}

CartesianGridGeometry = {
  'cell_dimensions': [800, 1, 1],
  'coordinates': {
    'lower': [-1.50, -0.015, -0.015],
    'upper': [-0.00, +0.015, +0.015],
  }
}

PatchHierarchy = {
 'max_number_of_levels': 1, 
 'refine_ratio': [2, 1, 1],
 'blocking_factor': [8, 1, 1],
 'max_grid_size': [120, 1, 1],
}

PressureValveBoundary = {
  'efficiency': 1.0,
  'open_at_interval': 0.04,
  'offset': 0,
  'fuel_measurement_position': -0.1,
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
    'required_massflow': 120.0 / 3600.0,
    'surface_area': math.pi * 0.015 * 0.015
  }
}

IgniteDetonation = {
  'interval': 0.06,
  'measurement_position': -0.1,
  'equivalence_ratio_criterium': 0.9,
  'position': -0.8,
}

IsentropicPressureBoundary = {
  'outer_pressure': 101325.0,
  'side': 1
}

Output = { 
  'outputs': [{
    'type': 'HDF5',
    'path': './IdealGasMix.h5',
    'intervals': [5e-5],
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  }]
}
