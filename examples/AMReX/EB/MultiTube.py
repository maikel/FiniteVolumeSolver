cfl = 0.95
final_time = 1.0
max_cycles = 10
max_number_of_levels = 1
plenum_n_cells = 64


output = {
    'Plotfile': { 
    'type': 'Plotfile',
    'directory': 'MultiTube/Test/',
    'intervals': [1e-4],
  }
}

ignite = {
  'interval': 0.03333333,
  'measurement_position': -0.2,
  'equivalence_ratio_criterium': 0.8,
  'position': 0.8
}

valves = dict()
for i in range(6):
  valves[i] = {
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

valve0 = valves[0]
valve1 = valves[1]
valve2 = valves[2]
valve3 = valves[3]
valve4 = valves[4]
valve5 = valves[5]


