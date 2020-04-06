cfl = 0.8
final_time = 1.0
max_cycles = 0
max_number_of_levels = 1
plenum_n_cells = 32


grid = {
  'n_cells': [64, 64, 64],
  'max_number_of_levels': 1
}

output = {
  'Plenum_x0': { 
    'type': 'HDF5',
    'filename': 'Plenum_x0.h5',
    'intervals': [1e-5],
  }
}

ignite = {
  'interval': 0.03333333,
  'measurement_position': 1.8,
  'equivalence_ratio_criterium': 0.8,
  'position': 0.8
}

valves = dict()
for i in range(5):
  valves[i] = {
    'efficiency': 1.0,
    'open_at_interval': 0.03333333,
    'offset': 0.04 + i * 6e-3,
    'fuel_measurement_position': 1.8,
    'fuel_measurement_criterium': 0.8,
    'pressure_value_which_opens_boundary': 101325.0,
    'pressure_value_which_closes_boundary': 3.0e5,
    'oxygen_measurement_position': 1.0,
    'oxygen_measurement_criterium': 0.1,
    'equivalence_ratio': 1.0
  }

valve0 = valves[0]
valve1 = valves[1]
valve2 = valves[2]
valve3 = valves[3]
valve4 = valves[4]


