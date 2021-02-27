import math

tube_n_cells = 256
# plenum_x_n_cells = 128
tube_blocking_factor = 8
plenum_blocking_factor = 8

massflow_cor = 0.24 #%MASSFLOW%
diffusion = 3.0 #%Diffusion%
outputPath = 'SEC_Plenum_Arrhenius' # '%OUTPUT%'

mode = 3 #%MODE%
boundary_condition = 'TurbineMassflowBoundaries' # '%BOUNDARY_CONDITION%'

n_level = 1

n_tubes = 6
r_tube = 0.015

D = 2.0 * r_tube

r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.2
r_tube_center = 2.0 * r_inner
alpha = 2.0 * math.pi / n_tubes

inlet_length = 6 * D # [m]

plenum_y_lower = - 0.5
plenum_y_upper = + 0.5
plenum_y_length = plenum_y_upper - plenum_y_lower

TVolRPlen = 20.0 * D
plenum_x_upper = TVolRPlen / plenum_y_length
plenum_x_lower = -0.3
plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_length = plenum_x_upper # [m]
tube_length = 1.0 # [m]

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

# tube_n_cells = plenum_x_n_cells * tube_over_plenum_length_ratio
# tube_n_cells -= tube_n_cells % tube_blocking_factor
# tube_n_cells = int(tube_n_cells)

# y0s = [-1.0/3.0, 0.0, +1.0/3.0]
y0s = [0.0]
n_tubes = len(y0s)

RunOptions = {
  'cfl': 0.1,# / float(tube_n_cells / 64),
  'final_time': 1.0,
  'max_cycles': -1,
  'do_backup': 0
}

LogOptions = {
  'file_template': 'Test-{rank}.txt',
  'channel_blacklist': ['TurbineMassflowBoundary']
}

FluxMethod = {
  # HLLEM, HLLEM_Larrouturou
  'base_method': 'HLLEM_Larrouturou',
  # Upwind, MinMod, VanLeer
  'limiter': 'MinMod',
  # Conservative, Primitive, Characteristics
  'reconstruction': 'Characteristics'
}

R = 1.0
gamma = 1.4

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
  'mul': diffusion
}
R_ref = 287.
p_ref = 101325.
T_ref = 300.
L_ref = 1.0
rho_ref = p_ref / T_ref / R_ref
u_ref = math.sqrt(p_ref / rho_ref)
t_ref = L_ref / u_ref
# ud->Msq =  u_ref*u_ref / (R_gas*T_ref);

p0 = 2.0
rho0 = math.pow(p0, 1.0 / gamma)
T0 = p0 / rho0
p = 0.95 * p0
T = T0 + ArrheniusKinetics['Q'] * (gamma - 1.0)
rho = p / T

# checkpoint = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Debug/Checkpoint/000000005'
checkpoint = ''


def Area(xi):
  A0  = 1.0
  A1  = 4.0 # Reference: 3.0; best: 4.0 
  xi0 = 0.0625 - 1.0
  xi1 = 0.75 - 1.0   # Reference: 0.5; best: 0.75
  Ax = 1.0 if xi < xi0 else A0 + (A1-A0)*(xi-xi0)/(xi1-xi0) if xi < xi1 else A1
  return Ax


ControlOptions = {
  'Q': ArrheniusKinetics['Q'],
  'rpmmin': (10000.0 / 60.) * t_ref,
  'rpmmax': (60000.0 / 60.) * t_ref,
  'initial_turbine_pressure': p,
  'initial_turbine_temperature': T,
  'target_pressure_compressor' : 6.0,
  'checkpoint': checkpoint,
  # Tube surface
  'surface_area_tube_inlet': (n_tubes * Area(-1.0) * D),
  'surface_area_tube_outlet': (n_tubes * Area(0.0) * D),
  # Turbine volumes and surfaces
  'volume_turbine_plenum': TVolRPlen,
  'surface_area_turbine_plenum_to_turbine': plenum_y_length,
  # Compressor volumes and surfaces
  'volume_compressor_plenum': TVolRPlen,
  'surface_area_compressor_to_compressor_plenum': 8.0 * D,
}

def ToCellIndex(x, xlo, xhi, ncells):
  xlen = xhi - xlo
  x_rel = (x - xlo) / xlen
  i = int(x_rel * ncells)
  return i

Plenum = {
  'checkpoint': checkpoint,
  'GridGeometry': {
    'cell_dimensions': [plenum_x_n_cells, plenum_y_n_cells, plenum_z_n_cells],
    'coordinates': {
      'lower': [plenum_x_lower, plenum_y_lower, plenum_z_lower],
      'upper': [plenum_x_upper, plenum_y_upper, plenum_z_upper],
    },
    # 'periodicity': [0, 1, 0]
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
    'scratch_gcw': 2,
    'flux_gcw': 0,
  },
  'FluxMethod': FluxMethod,
  'InletGeometries': [{
    'r_start': 4.0 * r_tube,
    'r_end': 4.0 * r_tube,
    'y_0': y_0,
    'height': inlet_length,
    'angle':  30. / 180. * math.pi
  } for y_0 in y0s],
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

Plenum[boundary_condition] = {
  'boundary_section': { 
    'lower': [plenum_x_n_cells, -4, 0], 
    'upper': [plenum_x_n_cells + 3, plenum_y_n_cells + 3, 0] 
    },
  'mode': mode,
  'coarse_average_mirror_box': {
    'lower': [plenum_x_n_cells - 1, 0, 0],
    'upper': [plenum_x_n_cells - 1, plenum_y_n_cells - 1, 0]
  },
  'relative_surface_area': ControlOptions['surface_area_turbine_plenum_to_turbine'], 
  'massflow_correlation': massflow_cor, # 0.06 * 4.0,
}

def TubeCenterPoint(x0, y0):
  return [x0, y0, 0.0]

def LowerX(x0, y0):
  center = TubeCenterPoint(x0, y0)
  center[1] -= r_tube
  center[2] -= r_tube
  return center

def UpperX(x0, y0):
  center = TubeCenterPoint(x0, y0)
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

def PlenumMirrorBox(y0):
  return BoxWhichContains(DomainAroundPoint(TubeCenterPoint(plenum_x_lower, y0), [0.0, -0.2], [-0.9 * plenum_x_lower, +0.2]))


Tube_FluxMethod = FluxMethod
Tube_FluxMethod['area_variation'] = Area

Tubes = [{
  'checkpoint': checkpoint if checkpoint == '' else '{}/Tube_{}'.format(checkpoint, i),
  # 'buffer': 0.06,
  'initially_filled_x': 0.1,
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
    'scratch_gcw': 4,
    'flux_gcw': 0,
  },
} for (i, y_0) in enumerate(y0s)]

mode_names = ['cellwise', 'average_mirror_state', 'average_ghost_state', 'average_massflow']

tube_intervals = 0.005
plenum_intervals = 0.02

# outputPath = mode_names[Plenum['TurbineMassflowBoundaries'][0]['mode']]


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
  {
    'type': 'Plotfiles',
    #'directory': '/group/ag_klima/SFB1029_C01/SEC_Plenum/{}/Plotfiles/'.format(mode_names[Plenum['TurbineMassflowBoundaries'][0]['mode']]),
    'directory': '{}/Plotfiles'.format(outputPath),
    'intervals': [0.005],
    #'frequencies': [1]
  }, {
    #'#type': 'Checkpoint',
    #'directory': 'SEC_Plenum_{}_{}/highres/Checkpoint/'.format(FluxMethod['base_method'], FluxMethod['limiter']),
    # 'intervals': [1.0],
    #'frequencies': [100]
  #}, {
   'type': 'CounterOutput',
  #  'intervals': [1/.0]
    'frequencies': [100]
  },
  {
    'type': 'Checkpoint',
    'intervals': [1.0],
    'directory': 'Checkpoint/'
  }
  ]
}
