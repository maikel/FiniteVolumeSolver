import math

# Controls the number of cells in the x-direction
nx = 1200
# The output files are written in this directory from the perspective of this 
# docker container
case = 'De1'
# This option defines the Mach number of the shock
# Default is 1.1
Mach_number = 1.61

############################################################################
# These variables are only auxiliary and are not directly used by the app 
#
base_path = './{}'.format(case)

def nCellsY(nCellsX, ratio, blocking_factor=8):
    base = int(nCellsX * ratio)
    ncells = base - base % blocking_factor
    return ncells

xlower = 0.0
xupper = 0.230000001
xlen = xupper - xlower

ylower = -0.05 + 0.005
yupper = +0.05 + 0.005
ylen = yupper - ylower
blocking_factor = 8
ny = nCellsY(nx, ylen / xlen, blocking_factor=blocking_factor)
n_level = 1
n_error_buf = 0 if n_level == 1 else 1

############################################################################
# These variables and dictionaries will be read by the application

#Equation = {
  # 'gamma': 1.4,
  # 'Rspec': 287
#}

RunOptions = {
  'cfl': 0.8,
  'final_time': 5e-4,
  'do_backup': False,
  'max_cycles': -1, # -1 means infinite and 0 means only initial condition
}

#checkpoint='/srv/public/Maikel/FiniteVolumeSolver/build_2D-RelWithDebugInfo/Divider2D/Checkpoint/000000174'

FluxMethod = {
  'reconstruction': 'Characteristics',
  'limiter':'VanLeer',
  'base_method': 'HLLEM_Larrouturou'
}

GridGeometry = {
  # This option tells the solver how many cells to use
  'cell_dimensions': [nx, ny, 1],
  # This option defines the coordinates of the computational domain
  'coordinates': {
    'lower': [xlower, ylower, +0.00],
    'upper': [xupper, yupper, +0.10],
  },
  # We have no periodicity in Divider
  'periodicity': [0, 0, 0]
}

# defines the number of ghost cells
scratch_gcw = 2 if n_level == 1 else 4

# defines the number of ghost faces
flux_gcw = scratch_gcw - 2

PatchHierarchy = {
  # Defines the number of refinement levels
  # 1 means no adaptive refinement
  'max_number_of_levels': n_level,
  # Set this to 1 if max_number_of_levles > 1
  'n_error_buf': [n_error_buf, n_error_buf, n_error_buf],
  # Defines the maximal size of a single patch
  # patches can by smaller thou
  'max_grid_size': [nx, ny, 1],
  # Do not change this, it defines how many ghost cells are needed for the
  # cut-cell geometry
  'ngrow_eb_level_set': max(scratch_gcw + 1, 9),
  # If a patch is completely covered by behind cut-cells remove them from the grid
  'remove_covered_grids': False,
  # Every patch size is a multiple of blocking_factor (in each direction)
  # i.e. blocking_factor = 8 means patches have the size of 8 or 16 or 24 or ...
  'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
  # How many cells are needed inbetween AMR levels
  'n_proper': 1,
}

# Adjust those paths to read the files that describe the wall boundaries for the Divider.
# The first and the last point need within each wall file need to be equal
wall_filenames = ['{}/wall_1.txt'.format(base_path),
                  '{}/wall_2.txt'.format(base_path),
                  '{}/wall_3.txt'.format(base_path)]

# Defines the 
Output = {
	'outputs': [
  # Write AMReX-output directories. This output preserves the AMR hierarchy 
  # It is readable from Paraview, VisIt and Python
  # {'type': 'Plotfiles', 'directory': 'Divider_DE5_{}/Plotfiles'.format(nx), 'intervals': [1e-5]},
  # Write simple HDF5 files
  # It is faster and writes a single grid without patches
  {'type': 'HDF5', 'path': 'Divider_{}_nx-{}_Ma-{}.h5'.format(case, nx, Mach_number), 'intervals': [1e-5]},
  # Print out timer statistics
  {'type': 'CounterOutput', 'intervals': [1e-4] },
  #{'type': 'Checkpoint', 'intervals': [1e-6] }
]
}
