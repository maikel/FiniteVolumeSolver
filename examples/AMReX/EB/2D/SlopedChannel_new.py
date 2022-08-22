import math

isRef = 0
factor = 1

final_time = 0.26

max_cycles = [10, 0]
mode = ['new', 'ref']
originC = [0.35, 0.35 + final_time]

RunOptions = {
  'cfl': 0.5,
  'final_time': final_time,
  'max_cycles': max_cycles[isRef] # -1 means infinite and 0 means only initial condition
}

def AlignForBlockingFactor(n, blocking_factor):
  return n - n % blocking_factor

blocking_factor = 1

n_cells_x = 50*factor
n_cells_y = AlignForBlockingFactor(35*factor, blocking_factor)
n_levels = 1

GridGeometry = {
  'cell_dimensions': [n_cells_x, n_cells_y, 1],
  'coordinates': {
    'lower': [+0.00, +0.00, +0.00],
    'upper': [+1.0, +0.7, +0.7],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_levels, 
  'n_error_buf': [1, 1, 1],
  'max_grid_size': [n_cells_x, n_cells_y, n_cells_x],
  'ngrow_eb_level_set': 9,
  'remove_covered_grids': False,
  'blocking_factor': [n_cells_x, n_cells_y, 1],
  'n_proper': 1,
  'hgrid_details': True
}

# limiter = "MinModLimiter"
limiter = "NoLimiter"
# limiter = "Upwind"
# reconstruction = "KBN"
reconstruction = "PrimitiveReconstruction"
# reconstruction = "ConservativeReconstruction"
# initial_function = "Sod"
initial_function = "Smooth"
# initial_function = "Linear"
# initial_function = "Constant"

rhoL = 1.0
uL = 1.0
pL = 1.0

rho0 = 1.0
u0 = 1.0
p0 = 1.0

# width = 0.3

rhoR = 1.0
uR = 1.0
pR = 1.0

origin = originC[isRef]
theta = math.pi * 30.0 / 180.0

Output = { 
  'outputs': [
    {
    'type': 'Plotfile',
    'directory': 'Debug_{}_{}x{}-{}-{}/Plot/'.format(reconstruction, n_cells_x, n_cells_y, limiter, initial_function),
    'intervals': [0.05],
    'frequencies': [1] 
  },
  # {
  #   'type': 'DebugOutput',
  #   'directory': 'Debug_{}_{}x{}-{}/'.format(reconstruction, n_cells_x, n_cells_y, limiter),
  #   'intervals': [1e-4],
  #   'frequencies': [1] 
  # },
  # {
  #   'type': 'HDF5',
  #   'path': 'ReferenceData/SlopedChannel_{}_{}_{}x{}-{}.h5'.format(mode[isRef], reconstruction, n_cells_x, n_cells_y, limiter),
  #   'intervals': [final_time],
  #   'frequencies': [1] 
  # }
  # {
    # 'type': 'CounterOutput',
    # 'intervals': [5.0e-5],
    #'frequencies': [1],
  # },
  ]
}
