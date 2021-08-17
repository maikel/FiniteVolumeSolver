import math

# This is the amount of cells in the x Direction for the plenum domain
# Given this and the total lengths will determine all cell sizes
plenum_x_nCells = 512

# This settings controls where and when the Riemann Problem 
# is applied in the Domain. 
# Further description is given below in ShockValveBoundary['schock_feedback']
shock_x_location = -0.05 
shock_mach_number =  2.62
shock_time = 100. #1.0e-02

# This controls the inflow velocity in the tube 
# The density in the domain is default 1.22 kg/m^3 
# so this correspond roughly to 48.0 m/s.
# Further description is given below in ShockValveBoundary['massflow_boundary']
r_tube = 0.016 # radius of the tube
required_massflow = 0.047 # [kg/s]
surface_area = math.pi * r_tube * r_tube # [m^2]

# The output files are written in this directory from the perspective of this 
# docker container
OutPut_BasePath = 'MB_Riemann_nx_{}_Ma_{}_mf_{}'.format(plenum_x_nCells, shock_mach_number, required_massflow)

# OutPut_BasePath = 'Div_Riemann_plenum_x_nCells_{}_Ma_{}_MF_{}'.format(plenum_x_nCells, shock_mach_number, required_massflow)


###############################################################
############      NUMERICAL SETTINGS      #####################
###############################################################

RunOptions = {
  # This is the CFL-condition. 0 < cfl < 1
  # Higher means larger time steps
  'cfl': 0.8,
  # The simulation stops upon reaching the final time
  'final_time': 2e-02, # 5e-1,
  'do_backup': False,
  # The simulation stops after doing max_cycles (many coarse time steps). 
  # Choose -1 for infinite many time steps.
  # Choose 0 for generating initial conditions.
  'max_cycles': -1,
  # Set the AxiSymmetric Source Term on or off.
  'AxiSymmetric' : True, # True / False
}

FluxMethod = {
  # TODO: add possible choices for all
  # Use characteristic variables at the reconstruction step
  'reconstruction': 'Characteristics',
  # Choose the family of slope limiter
  'limiter':'VanLeer',
  # Configure the underlying base method
  'base_method': 'HLLEM_Larrouturou'
}

# This may point to a directory that contains a checkpoint from a previous run
checkpoint = ''

###############################################################
#############      GEOMETRY SETTINGS      #####################
###############################################################
# These variables define the 1D Tube.
d_tube = 2.0 * r_tube

# This is the amount of the tube length that we will resolve in the 2D plenum domain
inlet_length = 3.0 * d_tube # [m]

tube_xlower = - ( 1.0 ) # [m]
tube_xupper = -inlet_length # [m]
# This is length of the computational domain containing the tube
tube_xlen = tube_xupper - tube_xlower #[m]

#-----------------------------------------------
# These variables define the 2D Plenum.

# Plenum starts logically at x = 0 and ends in plenum_x_upper
plenum_xlower = -inlet_length
plenum_xupper = 10.0 * d_tube
plenum_xlen = plenum_xupper - plenum_xlower

# Lower y-coordinate MUST be 0!!! Because of the Axi-Symmetric source term.
# Chose the upper coordinate such that you resolve everything you need.
plenum_ylower = 0.0
plenum_yupper = 10.0 * d_tube
plenum_ylen = plenum_yupper - plenum_ylower

# Adjust those paths to read the files that describe the wall boundaries for the Divider.
# The first and the last point need within each wall file need to be equal
wall_path = 'wall_files'
wall_filenames = [ '{}/wall.txt'.format(wall_path) ]

#-----------------------------------------------------

blocking_factor = 8
n_level = 1

def nCellsY(nCellsX, ratio, blocking_factor=8):
    base = int(nCellsX * ratio)
    ncells = base - base % blocking_factor
    return ncells

plenum_y_nCells = nCellsY(plenum_x_nCells, plenum_ylen / plenum_xlen, blocking_factor=blocking_factor)
tube_nCells = nCellsY(plenum_x_nCells, tube_xlen / plenum_xlen, blocking_factor=blocking_factor)

n_error_buf = 0 if n_level == 1 else 1

# defines the number of ghost cells
scratch_gcw = 4 if n_level == 1 else 4

# defines the number of ghost faces
flux_gcw = scratch_gcw - 2

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
    # This option tells the solver how many cells to use
    'cell_dimensions': [plenum_x_nCells, plenum_y_nCells, 1],
    # This option defines the coordinates of the computational domain
    'coordinates': {
      'lower': [plenum_xlower, plenum_ylower, 0.0],
      'upper': [plenum_xupper, plenum_yupper, 1.0],
    },
    # We have no periodicity
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    # Defines the number of refinement levels
    # 1 means no adaptive refinement
    'max_number_of_levels': n_level, 
    # Every patch size is a multiple of blocking_factor (in each direction)
    # i.e. blocking_factor = 8 means patches have the size of 8 or 16 or 24 or ...
    'blocking_factor': [blocking_factor, blocking_factor, blocking_factor],
    # Do not change this, it defines how many ghost cells are needed for the cut-cell geometry
    'ngrow_eb_level_set': max(scratch_gcw + 1, 9),
    # If a patch is completely covered by behind cut-cells remove them from the grid
    'remove_covered_grids': False,
    # How many cells are needed inbetween AMR levels
    'n_proper': 1,
    # Set this to 1 if max_number_of_levels > 1
    'n_error_buf': [n_error_buf, n_error_buf, n_error_buf],
  },
  'IntegratorContext': {
    # defines the number of ghost cells
    'scratch_gcw': scratch_gcw,
    'flux_gcw': flux_gcw
  }
}


Tube = {
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [tube_nCells, 1, 1],
    'coordinates': {
      'lower': LowerX(-tube_xlen),
      'upper': UpperX(-inlet_length),
    },
    'periodicity': [0, 0, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [blocking_factor, 1, 1],
    'max_grid_size': [tube_nCells, tube_nCells, tube_nCells],
    'refine_ratio': [2, 1, 1],
    'n_proper': 1,
    'n_error_buf': [4, 0, 0]
  },
  'IntegratorContext': {
    'scratch_gcw': 2,
    'flux_gcw': 0
  }
}

###############################################################
#############      PHYSICAL SETTINGS      #####################
###############################################################

# Define the adiabatic exponent and the specific as constant.
# Both default values are taken from dry air at standard conditions.
Equation = {
  'gamma': 1.4, # Default 1.4 [-]
  'Rspec': 287.058 # Default 287.058 [J/kg/K]
}

# This defines the constant initial state in the whole domain.
InitialCondition = {
    'density': 1.22, # Default 1.22 [kg/m^3]
    'u_velocity': 0.0, # Default 0.0 [m/s]
    'v_velocity': 0.0, # Default 0.0 [m/s]
    'w_velocity': 0.0, # Default 0.0 [m/s] # only for 3D currently not implemented !!!
    'pressure': 101325.0 # Default 101325.0 [Pa]
}

# The Variable shock_xid 
# is the cell number where this location lives.
# Further settings for the Riemann problem can be 
# set below in the dictionary ShockValveBoundary['schock_feedback'].
shock_xid = int(plenum_x_nCells * (shock_x_location-plenum_xlower) / plenum_xlen)

# This options influence the boundary condition on the left side of the domain.
ShockValveBoundary = {
  # If the valve is open we have a constant mass flow in the domain.
  'massflow_boundary': {
    # this is the box where we average the state and 
    # compute the neccessary velocity to set the mass flow
    'coarse_inner_box': { 
      'lower': [0, 0, 0], 
      'upper': [1, 0, 0] 
    },
    # This controls the inflow velocity in the tube 
    'required_massflow': required_massflow, # Default 0.0 [kg / s]
    'surface_area': surface_area, # Default 1.0 [m^2]
    # The surface area should be somehow connected to the Noozle Geometry.
  },
  # If the below given shock_time is reached 
  # we set up a Riemann Problem according to this data.
  # After this is done the valve is automatically set 
  # to a transmissive state and the mass flow is stopped.
  'schock_feedback': {
    # the Mach number of the Riemann Problem
    'shock_mach_number': shock_mach_number, # Default 1.0
    'shock_time': shock_time, # Default 10.0 [s]
    # In this box the state is averaged and defines 
    # the right state from the Riemann problem (post shock state).
    # After this is done we compute from this the 
    # left state from the Riemann problem.
    'average_post_shock_box' : {
      'lower': [shock_xid-1, 0, 0], 
      'upper': [shock_xid+1, plenum_y_nCells, 0] 
    },
  }
}

####################
# Don't edit, this is needed for constructing the shock_feedback class.
schock_feedback = ShockValveBoundary['schock_feedback']

# some Dicts must include each other 
Tube['ShockValveBoundary'] = ShockValveBoundary
Tube['Equation'] = Equation
Plenum['Equation'] = Equation
Tube['InitialCondition']=InitialCondition
Plenum['InitialCondition']=InitialCondition
Tube['FluxMethod']=FluxMethod
Plenum['FluxMethod']=FluxMethod
####################

###############################################################
##############      OUTPUT SETTINGS      ######################
###############################################################

LogOptions = {
  # this controls where the log file is written
  'file_template': '{}/0000.log'.format(OutPut_BasePath),
 }

# Defines the 
Output = {
	'outputs': [
  # Write AMReX-output directories. This output preserves the AMR hierarchy 
  # It is readable from Paraview, VisIt and Python
  {'type': 'Plotfiles',
    'directory': '{}/Plotfiles'.format(OutPut_BasePath),
    'intervals': [1e-4]
  },
  
  # Write simple HDF5 files
  # It is faster and writes a single grid without patches
  {'type': 'HDF5', 
    'path': '{}/Tube.h5'.format(OutPut_BasePath), 
    'intervals': [1e-4]
  },

  {'type': 'HDF5', 
    'path': '{}/Plenum.h5'.format(OutPut_BasePath), 
    'intervals': [1e-4]
  },
  
  # Print out timer statistics
  {'type': 'CounterOutput', 
    'intervals': [1e-2] 
  },
  
  # write out Checkpoints
  {'type': 'Checkpoint', 
    'directory': '{}/Checkpoint'.format(OutPut_BasePath),
    'intervals': [1e-3]  
  }
]
}
