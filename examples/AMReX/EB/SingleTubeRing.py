cfl = 0.8
final_time = 1.0
max_cycles = 0

x_length = 0.36
x_lower = -0.03

y_length = 0.3

x_range = [x_lower, x_lower + x_length]
y_range = [-0.5 * y_length, 0.5 * y_length]
z_range = [-0.5 * y_length, 0.5 * y_length]

def nCellsY(nCellsX, ratio):
  base = int(nCellsX * ratio)
  ncells = base - base % 8
  return ncells

nx = 96
ny = nCellsY(nx, y_length / x_length)
nz = ny

grid = {
  'n_cells': [nx, ny, nz],
  'max_number_of_levels': 2,
  'x_range': x_range,
  'y_range': y_range,
  'z_range': z_range
}

r_tube = 0.015

plenum = {
  'temperature': 300,
  'outlet_radius': 0.004,
  'length': 0.25,
  'jump': 2 * r_tube
}

output = {
  'intervals': [1e-3, 1e-5],
}

ignite = {
  'interval': 0.03333333,
  'measurement_position': 1.8,
  'equivalence_ratio_criterium': 0.8,
  'position': 0.8
}

valve = {
  'efficiency': 1.0,
  'open_at_interval': 0.03333333,
  'offset': 0.0,
  'fuel_measurement_position': 1.8,
  'fuel_measurement_criterium': 0.8,
  'pressure_value_which_opens_boundary': 1.2E5,
  'pressure_value_which_closes_boundary': 3.0e5,
  'oxygen_measurement_position': 1.0,
  'oxygen_measurement_criterium': 0.1,
  'equivalence_ratio': 1.0,
  'outer_pressure': 1.3E5
}
