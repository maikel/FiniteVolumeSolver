import math

RunOptions = {
  'cfl': 0.4,
  'final_time': 0.0002,
  'max_cycles': -1, # -1 means infinite and 0 means only initial condition
}

def nCellsY(nCellsX, ratio, blocking_factor=8):
    base = int(nCellsX * ratio)
    ncells = base - base % blocking_factor
    return ncells

Mach_number = 1.1

xlower = 0.0
xupper = 0.2300001
xlen = xupper - xlower

ylower = -0.05 + 0.005
yupper = +0.05 + 0.005
ylen = yupper - ylower
blocking_factor = 8
nx = 2400
ny = nCellsY(nx, ylen / xlen, blocking_factor=blocking_factor)

n_level = 1

GridGeometry = {
  'cell_dimensions': [nx, ny, 1],
  'coordinates': {
    'lower': [xlower, ylower, +0.00],
    'upper': [xupper, yupper, +0.10],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_level, 
  'n_error_buf': [0, 0, 0],
  'max_grid_size': [nx, ny, 1024],
  'ngrow_eb_level_set': 5,
  'remove_covered_grids': True,
  'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
  'n_proper': 1,
}


base_path = '/srv/public/Maikel/FiniteVolumeSolver/examples/AMReX/EB/2D/'

wall_filenames = ['{}/wall_1.txt'.format(base_path),
                  '{}/wall_2.txt'.format(base_path),
                  '{}/wall_3.txt'.format(base_path)]

# max_cycles = 1

Output = {
	'outputs': [
  # {'type': 'Checkpoint', 'directory': 'Divider_DE5/Checkpoint', 'intervals': [1e-4]},
  {'type': 'Plotfiles', 'directory': 'Divider_DE5/Plotfiles', 'intervals': [1e-4]},
  {'type': 'HDF5', 'path': 'Divider_c24.h5', 'intervals': [1e-5]}
]
}
