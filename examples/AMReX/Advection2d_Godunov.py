RunOptions = {
  'cfl': 1.0,# should be between in (0, 1]
  'final_time': 2.0,
  'max_cycles': -1, # -1 means infinite and 0 means only initial condition
}

GridGeometry = {
 'cell_dimensions': [200, 200],
 'coordinates': {
   'lower': [-1.0, -1.0],
   'upper': [+1.0, +1.0],
 },
 'periodicity': [1, 1]
}

PatchHierarchy = {
 'max_number_of_levels': 2, 
 'n_error_buf': [4, 4]
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'Advection_Godunov/',
    'intervals': [0.1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  }]
}
