import math

x_len = 1.0
n_cells = 256
tube_blocking_factor = 8

outputPath = 'sec_c_vali'

n_tubes = 1
r_tube = 0.015
D = 2.0 * r_tube

Rspec = 287.0
gamma = 1.4


RunOptions = {
  'cfl': 0.1125,# / float(tube_n_cells / 64),
  'final_time': 300.0,
  'max_cycles': -1,
  'do_backup': 0
}

LogOptions = {
  'file_template': '{}/000.log'.format(outputPath),
  'channel_blacklist': ['TurbineMassflowBoundary']
}

FluxMethod = {
  # HLLEM, HLLEM_Larrouturou
  'base_method': 'HLLEM_Larrouturou',
  # Upwind, MinMod, VanLeer
  'limiter': 'VanLeer',
  # Conservative, Primitive, Characteristics
  'reconstruction': 'Characteristics'
}

Equation = {
  'Rspec': Rspec,
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
R_ref = 287.4
p_ref = 10000.0
T_ref = 300.
L_ref = 1.0
rho_ref = p_ref / T_ref / R_ref
u_ref = math.sqrt(p_ref / rho_ref)
t_ref = L_ref / u_ref
Ma_sq =  u_ref*u_ref / (R_ref*T_ref)

p0 = 2.0
rho0 = math.pow(p0, 1.0 / gamma)
T0 = p0 / rho0 
p = 0.95 * p0
#T = T0 + ArrheniusKinetics['Q'] * (gamma - 1.0)
T = 11.290743302923245
rho = p / T


# checkpoint = '/home/zenkechr/FVS_develop/FiniteVolumeSolver/build_2D-Release/oneTube/Checkpoint/002025006'
checkpoint = ''

def Area(xi):
  A0  = 1.0
  A1  = 4.0 # Reference: 3.0; best: 4.0 
  xi0 = 0.0625
  xi1 = 0.75   # Reference: 0.5; best: 0.75
  Ax = 1.0 if xi < xi0 else A0 + (A1-A0)*(xi-xi0)/(xi1-xi0) if xi < xi1 else A1
  return Ax


Tube_FluxMethod = FluxMethod
Tube_FluxMethod['area_variation'] = Area

Tubes = {
  'checkpoint': checkpoint if checkpoint == '' else '{}/Tube_{}'.format(checkpoint),
  'initially_filled_x': 0.0,
  'FluxMethod': Tube_FluxMethod,
  'GridGeometry': {
    'cell_dimensions': [n_cells, 1, 1],
    'coordinates': {
    'lower': [-0.0, -0.015, -0.015],
    'upper': [+1.00, +0.015, +0.015],
    },
    'periodicity': [0, 0, 0]
  },
  'DiffusionSourceTerm': DiffusionSourceTerm,
  'ArrheniusKinetics': ArrheniusKinetics,
  'PatchHierarchy': {
    'max_number_of_levels': 1, 
    'max_grid_size': [n_cells, n_cells, n_cells],
    'blocking_factor': [tube_blocking_factor, 1, 1],
    'refine_ratio': [2, 1, 1],
    'n_proper': 1,
    'n_error_buf': [4, 0, 0]
  },
  'IntegratorContext': {
    'scratch_gcw': 6,
    'flux_gcw': 2,
  },
  'InflowOptionsSEC' : {
    'SEC_buffer': 0.06, # 0.06
    'SEC_tti': 1.2, # 1.2
    'SEC_timin': 0.1, # 0.1
  }
}

tube_x_lower = Tubes['GridGeometry']['coordinates']['lower'][0]
tube_x_upper = Tubes['GridGeometry']['coordinates']['upper'][0]

ControlOptions = {
  'Q': ArrheniusKinetics['Q'],
  'rpmmin': (10000.0 / 60.) * t_ref,
  'rpmmax': (60000.0 / 60.) * t_ref,
  'initial_turbine_pressure': p,
  'initial_turbine_temperature': T,
  'target_pressure_compressor' : 6.0,
  'checkpoint': checkpoint,
  # Tube surface
  'surface_area_tube_inlet': (Area(tube_x_lower) * D) / D,
  'surface_area_tube_outlet': (Area(tube_x_upper) * D) / D,
  # Turbine volumes and surfaces
  'volume_turbine_plenum': 20.0 * D / D,
  'surface_area_turbine_plenum_to_turbine': 4.0 * D / D,
  # Compressor volumes and surfaces
  'volume_compressor_plenum': 20.0 * D / D,
  'surface_area_compressor_to_compressor_plenum': (8.0 * D) / D,
}

ControlFeedback = {
  'mdot_correlation': 0.24
}

tube_intervals = 0.005

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
    'type': 'CounterOutput',
    # 'intervals': [1.0],
    'frequencies': [10000],
  },
  {
    'type': 'Checkpoint',
    'intervals': [1.0],
    'directory': '{}/Checkpoint/'.format(outputPath)
  }
  ]
}
