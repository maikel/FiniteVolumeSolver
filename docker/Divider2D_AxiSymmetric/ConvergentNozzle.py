import math

# Those vairables help with scaling the length by tube radius and tube diamenter
r_tube = 0.015   # Radius
D = 2.0 * r_tube # Diameter

# This is the amount of the tube length that we will resolve in the 2D plenum domain
inlet_length = 3.0 * D # [m]
# 2.083 is some number from the experiments. This is the effective length of the tube from the experiments
tube_x_lower = - (2.083 - 0.5) # [m]
tube_x_upper = -inlet_length # [m]
# Plenum starts logically at x = 0 and ends in plenum_x_upper
plenum_x_lower = -inlet_length # [m]
plenum_x_upper = 10.0 * D # [m]
# Lower y-coordinate MUST be 0!!! Because of the Axi-Symmetric source term.
# Chose the upper coordinate such that you resolve everything you need.
# If the boundary is problematic, this might be chosen big and then use adaptive mesh refinement
plenum_y_lower = 0.0
plenum_y_upper = 10.0 * D

# This is the length of the computation domain that conatains the plenum
plenum_x_length = plenum_x_upper - plenum_x_lower #[m]
# This is length of the computational domain containing the tube
tube_x_length = tube_x_upper - tube_x_lower #[m]


# This is the amount of cells in the x Direction for the plenum domain
# Given this and the total lengths will determine all cell sizes
plenum_x_n_cells = 256

# List of paths to ASCII files that contain points that form some kind of boundary
# Each file has the format
# <X1> <Y1>
# <X2> <Y2>
#  .     .
#  .     .
#  .     .
# <Xn> <Yn>
# <X1> <Y1>
# where <X1>, <Y1> are floating point numbers that describe the x- and y-coordinates of points.
# Each point in the sequence (Xi, Yi) will be connected with its neighbor points with a straight line.
# The last and first line need to be equal.
# The resulting boundary is the convex polygon spanned by those points.
# Each file in the list describes one body and the matehmatical union will be taken, i.e. intersections are possible.
wall_filenames = ['ConvergentNozzle/wall.txt']
# wall_filenames = ['Divergent_With_Body/wall.txt', 'Divergent_With_Body/body.txt']

# Number of refinement levels
n_level = 1

# This may point to a directory that contains a checkpoint from a previous run
checkpoint = ''

PressureValveBoundary = {
  # The offset says how long in the beginning until a change to fuel happens
  'change_to_fuel_time_offset': 10e-3, # [s]
  # This is the time interval between two cycles
  'change_to_fuel_at_interval': 30e-3, # [s]
  # This pressure value in Pascal [Pa] at which we change the inflow boudnary from fuel to air
  'pressure_value_which_closes_boundary': 3.0e5, # [Pa]
  # This controls the equivalence ratio of the fuel inflow condition
  'equivalence_ratio': 1.0,
  'massflow_boundary': {
    'coarse_inner_box': { 
      'lower': [0, 0, 0], 
      'upper': [1, 0, 0] 
    },
    'side': 0,
    'direction': 0,
    # This controls the inflow velocity in the tube 
    # This is taken from the PDC experiments and results in a flow velocity of roughly 40 m/s
    'required_massflow': 0.2 / 6.0, # [kg / s]
    'surface_area': math.pi * r_tube * r_tube
  }
}

fill_fraction = 0.7
IgniteDetonation = {
  # The minimal time duration between two ignitions
  'ignite_interval': 60e-3,
  # The x position at which a temperature hat will be put that ignites the fuel
  'ignite_position': -tube_x_length + 0.2,
  # The time offset for the first ignition. This ignition is disabled before this offset.
  'offset': 0.0,
  # Ignite the detonation if the equivalence ratio criterium is satisfied at measurement position.
  # This condition is only tested if the above ignition time constrictions are satisfied.
  # i.e. after the initial time offset and after we waited for the minimal duration between two ignitions.
  'measurement_position': -tube_x_length * (1.0 - fill_fraction),
  'equivalence_ratio_criterium': 0.9,
}

RunOptions = {
  # This is the CFL-condition. 0 < cfl < 1
  # Higher means larger time steps
  'cfl': 0.8,
  # The simulation stops upon reaching the final time
  'final_time': 0.04,
  # The simulation stops after doing max_cycles many coarse time steps. 
  # Choose -1 for infine many time steps.
  # Choose 0 for generating initial conditions.
  'max_cycles': -1
}

Output = { 
  'outputs': [
  # The next output saves the plenum domain data into a single HDF5 file
  # Change 'path' to the destination path of the file.
  # 'intervals' is the time interval at which the data is written.
  # The simulation will try to hit this interval even if small time steps are neccessary.
  {
    'type': 'HDF5',
    'which_block': 0,
    'path': 'ConvergentNozzle/Plenum.h5',
    'intervals': [1e-4]
  },
  # The next output saves the tube domain data into a single HDF5 file
  # Change 'path' to the destination path of the file.
  # 'intervals' is the time interval at which the data is written.
  # The simulation will try to hit this interval even if small time steps are neccessary.
  {
    'type': 'HDF5',
    'which_block': 1,
    'path': 'ConvergentNozzle/Tube.h5',
    'intervals': [1e-5]
  },
  # The next output saves the AMReX Plotfiles that are readable in Paraview or VisIt.
  # Change 'directory' to the destination directory. 
  # This output writes multiple directories called 'pltxxxxxx'. 
  # 'intervals' is the time interval at which the data is written.
  # The simulation will try to hit this interval even if small time steps are neccessary.
  {
    'type': 'Plotfiles',
    'directory': 'ConvergentNozzle/Plotfiles/',
    'intervals': [1e-4],
    # 'frequencies': [1]
  },
  # The next output saves the AMReX Plotfiles that are readable in Paraview or VisIt.
  # Change 'directory' to the destination directory. 
  # This output writes multiple directories called 'pltxxxxxx'. 
  # 'intervals' is the time interval at which the data is written.
  # The simulation will try to hit this interval even if small time steps are neccessary.
  {
    'type': 'Checkpoint',
    'directory': 'ConvergentNozzle/Checkpoint/',
    'intervals': [1e-3],
    # 'frequencies': [1]
  },
  # The next output writes performance counter output after frequencies many time steps
  {
    'type': 'CounterOutput',
    'frequencies': [200]
  }]
}

# All patch sizes are a multiple of the blocking factor
# i.e. plenum blocking factor = 8 means patches of size 8, 16, 24 etc...
tube_blocking_factor = 32
plenum_blocking_factor = 32

# The number of 
tube_over_plenum_length_ratio = tube_x_length / plenum_x_length
tube_n_cells = plenum_x_n_cells * tube_over_plenum_length_ratio
tube_n_cells -= tube_n_cells % tube_blocking_factor
tube_n_cells = int(tube_n_cells)

# This is the maximal patch size. If chosen small it generates alot of small patches.

plenum_y_length = plenum_y_upper - plenum_y_lower
plenum_y_over_x_ratio = plenum_y_length / plenum_x_length
plenum_y_n_cells = plenum_x_n_cells * plenum_y_over_x_ratio
plenum_y_n_cells -= plenum_y_n_cells % plenum_blocking_factor
plenum_y_n_cells = int(plenum_y_n_cells)

# Those functions help to compute certain coordinates regarding the combustion tube
def TubeCenterPoint(x0):
  return [x0, 0.0, 0.0]

def LowerX(x0):
  center = TubeCenterPoint(x0)
  return center

def UpperX(x0):
  center = TubeCenterPoint(x0)
  center[1] += r_tube
  center[2] += r_tube
  return center

Plenum = {
  'wall_filenames': wall_filenames,
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [plenum_x_n_cells, plenum_y_n_cells, 1],
    'coordinates': {
      'lower': [plenum_x_lower, plenum_y_lower, 0.0],
      'upper': [plenum_x_upper, plenum_y_upper, 1.0],
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [plenum_blocking_factor, plenum_blocking_factor, 1],
    'max_grid_size': [plenum_blocking_factor, plenum_blocking_factor, 1],
    'ngrow_eb_level_set': 5,
    'remove_covered_grids': False,
    'n_proper': 1,
    'n_error_buf': [0, 0, 0]
  },
  'IntegratorContext': {
    'scratch_gcw': 2,
    'flux_gcw': 0
  }
}


Tube = {
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [tube_n_cells, 1, 1],
    'coordinates': {
      'lower': LowerX(-tube_x_length),
      'upper': UpperX(-inlet_length),
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [tube_blocking_factor, 1, 1],
    'max_grid_size': [tube_n_cells, tube_n_cells, tube_n_cells],
    'refine_ratio': [2, 1, 1],
    'n_proper': 1,
    'n_error_buf': [4, 0, 0]
  },
  'PressureValveBoundary': PressureValveBoundary,
  'IntegratorContext': {
    'scratch_gcw': 2,
    'flux_gcw': 0
  }
}
