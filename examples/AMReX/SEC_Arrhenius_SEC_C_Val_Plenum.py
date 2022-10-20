import math

tube_length = 1.0
n_cells = 256 # 192 256 320 384
tube_blocking_factor = 8

n_tubes = 1
r_tube = 0.015
D = 2.0 * r_tube
y_0 = 0.0

TVolRPlen = 1.0 * 20.0 * D #20.0 * D

## InflowOptions for SEC Tube
## this values influence the SEC mode and are set under Tubes Dict
# fuel retardation time (controls how long the air buffer is)
SEC_buffer = 0.06 # default value 0.06
# maximum ignition delay time for the fuel
SEC_tti = 1.2 # default value 1.2
# minimum ignition delay time for the fuel
SEC_timin = 0.1 # default value 0.1

# parameter for the diffusorpart from the inflow tube
diffusorStart = 1.0 * 0.25 # start point x-axis
diffusorEnd = 0.75 # end point x-axis
offset = -tube_length # tube begins at x=-1.0
A0=1.0 # surface Area befor diffusor
A1=4.0 # surface Area after diffusor

outputPath = 'sec_c_vol{}_tx{}_xi0_{}_buf{}'.format(TVolRPlen/D, n_cells, diffusorStart, SEC_buffer)

Rspec = 287.0 # why real value??? in EB Case we use R=1.0
gamma = 1.4

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
  xi0 = diffusorStart
  xi1 = diffusorEnd
  
  xi0 += offset
  xi1 += offset # Reference: 0.5; best: 0.75
  Ax = 1.0 if xi < xi0 else A0 + (A1-A0)*(xi-xi0)/(xi1-xi0) if xi < xi1 else A1
  return Ax

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

Tube_FluxMethod = FluxMethod
Tube_FluxMethod['area_variation'] = Area

Tubes = {
  'checkpoint': checkpoint if checkpoint == '' else '{}/Tube_{}'.format(checkpoint),
  'initially_filled_x': 0.4,
  'FluxMethod': Tube_FluxMethod,
  'GridGeometry': {
    'cell_dimensions': [n_cells, 1, 1],
    'coordinates': {
    'lower': LowerX(-tube_length, y_0),
    'upper': UpperX(0.0, y_0),
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
    'SEC_buffer': SEC_buffer,
    'SEC_tti': SEC_tti,
    'SEC_timin': SEC_timin,
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
  'volume_turbine_plenum': TVolRPlen / D,
  'surface_area_turbine_plenum_to_turbine': 4.0 * D / D,
  # Compressor volumes and surfaces
  'volume_compressor_plenum': TVolRPlen / D,
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
    'intervals': [10.0],
    'directory': '{}/Checkpoint/'.format(outputPath)
  }
  ]
}
