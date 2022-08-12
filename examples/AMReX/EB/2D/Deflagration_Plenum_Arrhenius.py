import math

tube_n_cells = 192 # 192 256 320 384
# plenum_x_n_cells = 128
tube_blocking_factor = 8
plenum_blocking_factor = 32 #8

mode = 3 #%MODE%
boundary_condition = 'TurbineMassflowBoundaries' # '%BOUNDARY_CONDITION%'

n_level = 1

# tube_y0s = [-1.0/3.0, 0.0, +1.0/3.0]
tube_y0s = [0.0]
tube_length = 1.0 # [m]

n_tubes = len(tube_y0s)
r_tube = 0.015
D = 2.0 * r_tube

# parameter for the diffusorpart from the inflow tube
diffusorStart = 1.0 * 0.25 # start point x-axis
diffusorEnd = 0.75 # end point x-axis
offset=-tube_length # tube begins at x=-1.0
A0=1.0 # surface Area befor diffusor
A1=4.0 # surface Area after diffusor

# calculate plenum geometry
magic_z_length = 1.0 # should be replaced when switch to 3d!!!

# normally 3.0 * D # [m] # in old slanted case 10.0D
inlet_length = 3.0 * D 

plenum_y_lower = - 1.2
plenum_y_upper = + 1.2
plenum_y_length = plenum_y_upper - plenum_y_lower

TVolRPlen = 1.0 * 20.0 * D #20.0 * D
plenum_x_upper = TVolRPlen / plenum_y_length / magic_z_length
plenum_x_lower = -inlet_length
plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_length = plenum_x_upper - 0.0 # [m]

plenum_max_grid_size = max(plenum_blocking_factor, 1024)

plenum_domain_length = plenum_length + inlet_length
tube_domain_length = tube_length - inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length
plenum_over_tube_length_ratio = 1.0 / tube_over_plenum_length_ratio

plenum_x_n_cells = tube_n_cells * plenum_over_tube_length_ratio
plenum_x_n_cells -= plenum_x_n_cells % plenum_blocking_factor
plenum_x_n_cells = int(plenum_x_n_cells)

plenum_z_lower = plenum_y_lower
plenum_z_upper = plenum_y_upper
plenum_z_length = plenum_z_upper - plenum_z_lower

plenum_y_over_x_ratio = plenum_y_length / plenum_x_length

plenum_y_n_cells = plenum_x_n_cells * plenum_y_over_x_ratio
plenum_y_n_cells -= plenum_y_n_cells % plenum_blocking_factor
plenum_y_n_cells = int(plenum_y_n_cells)

plenum_z_over_x_ratio = plenum_z_length / plenum_x_length

plenum_z_n_cells = plenum_x_n_cells * plenum_z_over_x_ratio
plenum_z_n_cells -= plenum_z_n_cells % plenum_blocking_factor
plenum_z_n_cells = int(plenum_z_n_cells)


tube_dx = tube_domain_length/ tube_n_cells
plenum_dx = plenum_domain_length / plenum_x_n_cells
plenum_dy = plenum_y_length / plenum_y_n_cells
plenum_dz = plenum_z_length / plenum_z_n_cells

# formatter = (4*"{:>12}")

# print(formatter.format('', 'ncells', 'length', 'delta_x'))
# print(formatter.format('Tube_x', tube_n_cells, tube_domain_length, round(tube_dx, 5)))
# print(formatter.format('plenum_x', plenum_x_n_cells, plenum_domain_length, round(plenum_dx, 5)))
# print(formatter.format('plenum_y', plenum_y_n_cells, plenum_y_length, round(plenum_dy, 5)))
# print(formatter.format('plenum_z', plenum_z_n_cells, plenum_z_length, round(plenum_dz, 5)))


outputPath = 'defl_vol{}_y{}_tx{}_px{}_py{}_xi0_{}'.format(TVolRPlen/D, plenum_y_upper, tube_n_cells, plenum_x_n_cells, plenum_y_n_cells, diffusorStart)

RunOptions = {
  'cfl': 0.1125,# / float(tube_n_cells / 64),
  'final_time': 400.0,
  'max_cycles': -1,
  'do_backup': 0
}

LogOptions = {
  'file_template': '{}/000.log'.format(outputPath),
  'channel_blacklist': ['TurbineMassflowBoundary']
}

InputFileOptions = {
  'copy_input_file' : 1,
  'file_template': '{}/{}'.format(outputPath, 'inputfile.py'),
}

FluxMethod = {
  # HLLEM, HLLEM_Larrouturou
  'base_method': 'HLLEM_Larrouturou',
  # Upwind, MinMod, VanLeer
  'limiter': 'VanLeer',
  # Conservative, Primitive, Characteristics
  'reconstruction': 'Characteristics'
}

R = 1.0 # non dimensionalized specific gas constant
gamma = 1.4 # adiabitic exponent (for air)

Equation = {
  'Rspec': R,
  'gamma': gamma
}

T1ovT0 = 1.889
EaovRT1  = 10.0

ArrheniusKinetics = {
  'Q': 10.0/(gamma-1.0),
  'EA': EaovRT1 * T1ovT0,
  'B': (0.11/1.25)*math.exp(EaovRT1*(1.0-T1ovT0)),
  'T_switch': 1.6,
}

DiffusionSourceTerm = {
  'mul': 3.0
}

#-----------------------------------------------
# parameters for non dimensionalizing
R_ref = 287.4 # [J/kg/K] should be dry air
p_ref = 1.0e5 # [Pa]
T_ref = 300. # [K]
L_ref = 1.0 # [m]
rho_ref = p_ref / T_ref / R_ref
u_ref = math.sqrt(p_ref / rho_ref)
t_ref = L_ref / u_ref
# ud->Msq =  u_ref*u_ref / (R_gas*T_ref);

#-----------------------------------------------
# initial parameters for Plenum (used in Plenum dictionary below)
p0 = 2.0
rho0 = math.pow(p0, 1.0 / gamma)
T0 = p0 / rho0

p = 0.95 * p0
# T = T0 + ArrheniusKinetics['Q'] * (gamma - 1.0)
T = 11.290743302923245 # this value is used in Klein's Code
rho = p / T

#-----------------------------------------------
# checkpoint = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Debug/Checkpoint/000000005'
checkpoint = ''

def Area(xi):
  xi0 = diffusorStart
  xi1 = diffusorEnd
  
  xi0 += offset
  xi1 += offset # Reference: 0.5; best: 0.75
  Ax = 1.0 if xi < xi0 else A0 + (A1-A0)*(xi-xi0)/(xi1-xi0) if xi < xi1 else A1
  return Ax

# should be Tubes['GridGeometry']['coordinates']['lower'][0]
surface_area_SingleTube_inlet = Area(-1.0)
# should be Tubes['GridGeometry']['coordinates']['upper'][0]
surface_area_SingleTube_outlet = Area(0.0)

ControlOptions = {
  'Q': ArrheniusKinetics['Q'],
  'rpmmin': (10000.0 / 60.) * t_ref,
  'rpmmax': (60000.0 / 60.) * t_ref,
  'initial_turbine_pressure': p,
  'initial_turbine_temperature': T,
  'target_pressure_compressor' : 6.0,
  'checkpoint': checkpoint,
  # Tube surface
  'surface_area_tube_inlet': (n_tubes * surface_area_SingleTube_inlet * D) / D,
  'surface_area_tube_outlet': (n_tubes * surface_area_SingleTube_outlet * D) / D,
  # Turbine volumes and surfaces
  'volume_turbine_plenum': TVolRPlen / D,
  'surface_area_turbine_plenum_to_turbine': 4.0 * D / D,
  # Compressor volumes and surfaces
  'volume_compressor_plenum': TVolRPlen / D,
  'surface_area_compressor_to_compressor_plenum': (8.0 * D) / D,
}

def ToCellIndex(x, xlo, xhi, ncells):
  xlen = xhi - xlo
  x_rel = (x - xlo) / xlen
  i = int(x_rel * ncells)
  return i

def GetCenterPoint(x0, y0):
  return [x0, y0, 0.0]

def LowerX(x0, y0):
  center = GetCenterPoint(x0, y0)
  center[1] -= r_tube
  center[2] -= r_tube
  return center

def UpperX(x0, y0):
  center = GetCenterPoint(x0, y0)
  center[1] += r_tube
  center[2] += r_tube
  return center

def DomainAroundPoint(x0, lo, upper):
  xlo = [x0[0] + lo[0], x0[1] + lo[1]]
  xhi = [x0[0] + upper[0], x0[1] + upper[1]]
  return [xlo, xhi]

def BoxWhichContains(real_box):
  i0 = ToCellIndex(real_box[0][0], plenum_x_lower, plenum_x_upper, plenum_x_n_cells)
  iEnd = ToCellIndex(real_box[1][0], plenum_x_lower, plenum_x_upper, plenum_x_n_cells)
  j0 = ToCellIndex(real_box[0][1], plenum_y_lower, plenum_y_upper, plenum_y_n_cells)
  jEnd = ToCellIndex(real_box[1][1], plenum_y_lower, plenum_y_upper, plenum_y_n_cells)
  return { 'lower': [i0, j0, 0], 'upper': [iEnd, jEnd, 0] }

def BoxWhichContains_withGhostcells(real_box, gcw_x, gcw_y):
  i0 = ToCellIndex(real_box[0][0], plenum_x_lower, plenum_x_upper, plenum_x_n_cells)
  iEnd = ToCellIndex(real_box[1][0], plenum_x_lower, plenum_x_upper, plenum_x_n_cells)
  j0 = ToCellIndex(real_box[0][1], plenum_y_lower, plenum_y_upper, plenum_y_n_cells)
  jEnd = ToCellIndex(real_box[1][1], plenum_y_lower, plenum_y_upper, plenum_y_n_cells)
  return { 'lower': [i0-gcw_x[0], j0-gcw_y[0], 0], 'upper': [iEnd+gcw_x[1]-1, jEnd+gcw_y[1]-1, 0] }


Plenum = {
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [plenum_x_n_cells, plenum_y_n_cells, plenum_z_n_cells],
    'coordinates': {
      'lower': [plenum_x_lower, plenum_y_lower, plenum_z_lower],
      'upper': [plenum_x_upper, plenum_y_upper, plenum_z_upper],
    },
    'periodicity': [0, 1, 0]
  },
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'blocking_factor': [plenum_blocking_factor, plenum_blocking_factor, plenum_blocking_factor],
    'max_grid_size': [plenum_max_grid_size, plenum_max_grid_size, plenum_max_grid_size],
    'ngrow_eb_level_set': 9,
    'remove_covered_grids': False,
    'n_proper': 1,
    'n_error_buf': [0, 0, 0]
  },
  'IntegratorContext': {
    'scratch_gcw': 4,
    'flux_gcw': 2,
  },
  'FluxMethod': FluxMethod,
  'InletGeometries': [{
    'r_start': 4.0 * r_tube,
    'r_end': 4.0 * r_tube,
    'y_0': y_0,
  } for y_0 in tube_y0s],
  'InitialCondition': {
    'left': {
      'density': rho,
      'temperature': T,
      'pressure': p
    },
    'right': {
      'density': rho,
      'temperature': T,
      'pressure': p
    },
  }
}


Plenum_scratch_gcw = Plenum['IntegratorContext']['scratch_gcw']

# get Boundarybox for the Plenum (upper side in x-Direction)
plenum_y_midpoint = (plenum_y_upper + plenum_y_lower)/2.0
def PlenumBoundaryBox(y0):
  return BoxWhichContains_withGhostcells(DomainAroundPoint( GetCenterPoint(plenum_x_upper, y0), [0.0, -2.*D], [0.0, +2.*D]), 
            [0,Plenum_scratch_gcw], [0, 0])

## if plenum boundary goes over the whole y-range
# def PlenumBoundaryBox(y0):
#   return BoxWhichContains_withGhostcells(DomainAroundPoint( GetCenterPoint(plenum_x_upper, y0), [0.0, plenum_y_lower], [0.0, plenum_y_upper]), 
#             [0,Plenum_scratch_gcw], [Plenum_scratch_gcw, Plenum_scratch_gcw])


# get the Box for the Plenum where the massflow is accumulated (upper side in x-Direction)
def PlenumCoarseAverageMirrorBox(y0):
  return BoxWhichContains_withGhostcells(DomainAroundPoint( GetCenterPoint(plenum_x_upper, y0), [0.0, -2.*D], [0.0, +2.*D]), 
            [1,0], [0, 0])

Plenum[boundary_condition] = {
  # 'boundary_section': { 
  #   'lower': [plenum_x_n_cells, - Plenum_scratch_gcw, 0], 
  #   'upper': [plenum_x_n_cells + Plenum_scratch_gcw - 1, plenum_y_n_cells + Plenum_scratch_gcw - 1, 0] 
  #   },
  'boundary_section': PlenumBoundaryBox(plenum_y_midpoint),
  'mode': mode,
  # 'coarse_average_mirror_box': {
  #   'lower': [plenum_x_n_cells - 1, 0, 0],
  #   'upper': [plenum_x_n_cells - 1, plenum_y_n_cells - 1, 0]
  # },
  'coarse_average_mirror_box': PlenumCoarseAverageMirrorBox(plenum_y_midpoint),
  'relative_surface_area': ControlOptions['surface_area_turbine_plenum_to_turbine'], 
  'massflow_correlation': 0.06 * 4.0,
}

# get Mirrorbox for each Tube where they are connected with the plenum
def PlenumMirrorBox(y0):
  return BoxWhichContains(DomainAroundPoint(GetCenterPoint(plenum_x_lower, y0), [0.0, -2.0 * D], [inlet_length, +2.0 * D]))

# print(PlenumMirrorBox(tube_y0s[0]))
Tube_FluxMethod = FluxMethod
Tube_FluxMethod['area_variation'] = Area

Tubes = [{
  'checkpoint': checkpoint if checkpoint == '' else '{}/Tube_{}'.format(checkpoint, i),
  'initially_filled_x': 0.4,
  'FluxMethod': Tube_FluxMethod,
  'plenum_mirror_box': PlenumMirrorBox(y_0),
  'GridGeometry': {
    'cell_dimensions': [tube_n_cells, 1, 1],
    'coordinates': {
      'lower': LowerX(-tube_length, y_0),
      'upper': UpperX(-inlet_length, y_0),
    },
    'periodicity': [0, 0, 0]
  },
  'DiffusionSourceTerm': DiffusionSourceTerm,
  'ArrheniusKinetics': ArrheniusKinetics,
  'PatchHierarchy': {
    'max_number_of_levels': n_level, 
    'max_grid_size': [tube_n_cells, tube_n_cells, tube_n_cells],
    'blocking_factor': [tube_blocking_factor, 1, 1],
    'refine_ratio': [2, 1, 1],
    'n_proper': 1,
    'n_error_buf': [4, 0, 0]
  },
  'IntegratorContext': {
    'scratch_gcw': 6,
    'flux_gcw': 2,
  },
} for (i, y_0) in enumerate(tube_y0s)]

tube_intervals = 0.005
plenum_intervals = 0.02

Output = { 
  'outputs': [
    {
      'type': 'ControlOutput',
      'path': '{}/ControlState.h5'.format(outputPath),
      'intervals': [0.001],
      # 'frequencies': [1]
    },
  {
    'type': 'HDF5',
    'path': '{}/Tube0.h5'.format(outputPath),
    'which_block': 1,
    'intervals': [tube_intervals],
    # 'frequencies': [1]
  },
  {
    'type': 'HDF5',
    'path': '{}/Plenum.h5'.format(outputPath),
    'which_block': 0,
    'intervals': [plenum_intervals],
    # 'frequencies': [1]
  },
#  {
#    'type': 'Plotfiles',
#    'directory': '{}/Plotfiles/'.format(outputPath),
#    'intervals': [0.01],
#  },
  {
   'type': 'CounterOutput',
  #  'intervals': [1/.0]
    'frequencies': [10000]
  },
  {
    'type': 'Checkpoint',
    'intervals': [10.0],
    'directory': '{}/Checkpoint/'.format(outputPath)
  }
  ]
}
