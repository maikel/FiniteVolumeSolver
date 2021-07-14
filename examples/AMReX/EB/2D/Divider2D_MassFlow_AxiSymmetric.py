import math

# Controls the number of cells in the x-direction
nx = 128

# This settings controls where and when the Riemann Problem 
# is applied in the Domain. 
# Further description is given below in ShockValveBoundary['schock_feedback']
shock_x_location = -0.05 
shock_mach_number = 1.2
shock_time = 3.3e-03

# This controls the inflow velocity in the tube 
# The density in the domain is default 1.22 kg/m^3 
# so this correspond roughly to 35.0 m/s
# Further description is given below in ShockValveBoundary['massflow_boundary']
required_massflow = 35.0
surface_area = 1.0/1.22

# The output files are written in this directory from the perspective of this 
# docker container
OutPut_BasePath = 'Div_Riemann_nx_{}_Ma_{}'.format(nx, shock_mach_number)

###############################################################
############      NUMERICAL SETTINGS      #####################
###############################################################

RunOptions = {
  'cfl': 0.8,
  'final_time': 5e-3,
  'do_backup': False,
  'max_cycles': -1, # -1 means infinite and 0 means only initial condition
  # Set the AxiSymmetric Source Term on or off.
  'AxiSymmetric' : True, # True / False
}

FluxMethod = {
  'reconstruction': 'Characteristics',
  'limiter':'VanLeer',
  'base_method': 'HLLEM_Larrouturou'
}

###############################################################
#############      GEOMETRY SETTINGS      #####################
###############################################################
# These variables define the 2D Plenum.

xlower = -0.1
xupper = 0.1
xlen = xupper - xlower

# Lower y-coordinate MUST be 0!!! Because of the Axi-Symmetric source term.
# Chose the upper coordinate such that you resolve everything you need.
ylower = 0.0
yupper = 0.1

# Adjust those paths to read the files that describe the wall boundaries for the Divider.
# The first and the last point need within each wall file need to be equal
wall_path = '../examples/AMReX/EB/2D/Divider2D_MassFlow_AxiSymmetric/'
wall_filenames = ['{}/wall.txt'.format(wall_path),
                  # '{}/wall_2.txt'.format(wall_path),
                  ]

# These variables are only auxiliary and are not directly used by the app 
#
def nCellsY(nCellsX, ratio, blocking_factor=8):
    base = int(nCellsX * ratio)
    ncells = base - base % blocking_factor
    return ncells

ylen = yupper - ylower
blocking_factor = 8
ny = nCellsY(nx, ylen / xlen, blocking_factor=blocking_factor)
n_level = 1
n_error_buf = 0 if n_level == 1 else 1

GridGeometry = {
  # This option tells the solver how many cells to use
  'cell_dimensions': [nx, ny, 1],
  # This option defines the coordinates of the computational domain
  'coordinates': {
    'lower': [xlower, ylower, +0.00],
    'upper': [xupper, yupper, +0.10],
  },
  # We have no periodicity
  'periodicity': [0, 0, 0]
}

# defines the number of ghost cells
scratch_gcw = 4 if n_level == 1 else 4

# defines the number of ghost faces
flux_gcw = scratch_gcw - 2

PatchHierarchy = {
  # Defines the number of refinement levels
  # 1 means no adaptive refinement
  'max_number_of_levels': n_level,
  # Set this to 1 if max_number_of_levels > 1
  'n_error_buf': [n_error_buf, n_error_buf, n_error_buf],
  # Defines the maximal size of a single patch
  # patches can be smaller than these
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
shock_xid = int(nx * (shock_x_location-xlower) / xlen)

# This options influence the boundary condition on the left side of the domain.
ShockValveBoundary = {
  # If the valve is open we have a constant mass flow in the domain.
  'massflow_boundary': {
    # this is the box where we average the state and 
    # compute the neccessary velocity to set the mass flow
    'coarse_inner_box': { 
      'lower': [0, 0, 0], 
      'upper': [2, ny, 0] 
    },
    # This controls the inflow velocity in the tube 
    # The density in the domain is 1.22 kg/m^3 so this correspond roughly to 35.0 m/s
    'required_massflow': required_massflow, # Default 0.0 [kg / s]
    'surface_area': surface_area # Default 0.0 [m^2]
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
      'upper': [shock_xid+1, ny, 0] 
    },
  }
}

# Don't edit, this is needed for constructing the shock_feedback class.
schock_feedback = ShockValveBoundary['schock_feedback']

###############################################################
##############      OUTPUT SETTINGS      ######################
###############################################################

LogOptions = {
  'file_template': '{}/0000.log'.format(OutPut_BasePath),
 }

# Defines the 
Output = {
	'outputs': [
  # Write AMReX-output directories. This output preserves the AMR hierarchy 
  # It is readable from Paraview, VisIt and Python
  # {'type': 'Plotfiles', 'directory': '{}/Plotfiles'.format(OutPut_BasePath), 'intervals': [1e-5]},
  
  # Write simple HDF5 files
  # It is faster and writes a single grid without patches
  {'type': 'HDF5', 'path': '{}/Plenum.h5'.format(OutPut_BasePath, nx, shock_mach_number), 'intervals': [1e-5]},
  
  # Print out timer statistics
  {'type': 'CounterOutput', 'intervals': [1e-4] },
  
  # write out Checkpoints
  {'type': 'Checkpoint', 'directory': '{}/Checkpoint'.format(OutPut_BasePath),'intervals': [1e-3]  }
]
}
