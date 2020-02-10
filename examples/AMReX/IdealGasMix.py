RunOptions = {
  'cfl': 0.5,
  'final_time': 0.001,
}

CartesianGridGeometry = {
  'cell_dimensions': [15000, 1],
  'coordinates': {
    'lower': [-1.50, -0.015],
    'upper': [-0.00, +0.015],
  },
  'periodicity': [0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': 1, 
  'refine_ratio': [2, 1],
  'blocking_factor': [8, 1],
  'max_grid_size': [120, 1],
  'periodicity': [0, 0]
}

TemperatureRamp = {
  'left': {
    'pressure': 1.0 * 101325.0
  }, 'right': {
    'pressure': 1.0 * 101325.0
  },
  'fill_fraction': -0.3,
  'ramp_width': 0.01
}

IsentropicPressureBoundary = {
  'outer_pressure': 101325.0,
  'side': 1
}

Output = { 
  'outputs': [{
    'type': 'HDF5',
    'path': './IdealGasMix.h5',
    'intervals': [1e-5],
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  }]
}
