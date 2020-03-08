n_level = 1

tube_length = 0.15 # [m]
tube_radius = 0.015 # [m]
plenum_radius= 0.15 # [m]

blocking_factor = 8
max_grid_size = max(blocking_factor, 64)

eps = 0.0 # 3e-4 + 1e-5 

x_upper = plenum_radius + 0.01
x_lower = -1.0
x_length = x_upper - x_lower

y_upper = 0.02 + eps
y_lower = -2*plenum_radius - 0.01 + eps
y_length = y_upper - y_lower

z_upper = 0.0 # plenum_radius + 0.01 + eps
z_lower = -plenum_radius - 0.01 + eps
z_length = z_upper - z_lower

n_cells_y = 128
n_cells_y -= n_cells_y % blocking_factor
n_cells_y = int(n_cells_y)

x_over_y_ratio = x_length / y_length
z_over_y_ratio = z_length / y_length

n_cells_z = n_cells_y * z_over_y_ratio
n_cells_z -= n_cells_z % blocking_factor
n_cells_z = int(n_cells_z)

n_cells_x = n_cells_y * x_over_y_ratio
n_cells_x -= n_cells_x % blocking_factor
n_cells_x = int(n_cells_x)

Mach_number = 1.9

RunOptions = {
  'cfl': 0.5,
  'final_time': 0.004,
  'max_cycles': -1,
  'do_backup': 0
}

GridGeometry = {
  'cell_dimensions': [n_cells_x, n_cells_y, n_cells_z],
  'coordinates': {
    'lower': [x_lower, y_lower, z_lower],
    'upper': [x_upper, y_upper, z_upper],
  }
}

InitialCondition = {
  'left': {
    'moles': 'H2O:0.2, N2:0.8',
    'temperature': 1800,
    'pressure': 4.0 * 101325.0
  }, 'right': {
    'moles': 'O2:0.2, N2:0.8',
    'temperature': 300,
    'pressure': 101325.0
  }
}

PatchHierarchy = {
  'max_number_of_levels': n_level, 
  'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
  'max_grid_size': [max_grid_size, max_grid_size, max_grid_size],
  'ngrow_eb_level_set': 5,
  'remove_covered_grids': True
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'Kugelplenum/',
    'intervals': [1e-5],
    # 'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [25]
  }]
}
