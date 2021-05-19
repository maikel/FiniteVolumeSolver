import math

tube_n_cells = 128
# plenum_x_n_cells = 128
tube_blocking_factor = 8
plenum_blocking_factor = 8

angle_degree = 0.
angle = angle_degree / 180. * math.pi
outputPath = 'slanted_Plenum_test_angle{}_vol40_y0.48'.format(int(angle_degree))

mode = 3 #%MODE%
boundary_condition = 'TurbineMassflowBoundaries' # '%BOUNDARY_CONDITION%'

n_level = 1

y0s = [0.0]
n_tubes = len(y0s)

r_tube = 0.015
D = 2.0 * r_tube

inlet_length = 10. * D # [m]
tube_length = 1.0 # [m]
tube_domain_length = tube_length - inlet_length

plenum_max_grid_size = max(plenum_blocking_factor, 1024)
plenum_y_lower = - 0.48
plenum_y_upper = + 0.48
plenum_y_length = plenum_y_upper - plenum_y_lower

TVolRPlen = 2.0 * 20.0 * D
plenum_x_upper = TVolRPlen / plenum_y_length

plenum_length = plenum_x_upper # [m]
plenum_domain_length = plenum_length + inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length
plenum_over_tube_length_ratio = 1.0 / tube_over_plenum_length_ratio

plenum_x_n_cells = tube_n_cells * plenum_over_tube_length_ratio
plenum_x_n_cells -= plenum_x_n_cells % plenum_blocking_factor
plenum_x_n_cells = int(plenum_x_n_cells)

plenum_dx = plenum_domain_length / plenum_x_n_cells
plenum_x_lower = -inlet_length #-plenum_dx
plenum_x_length = plenum_x_upper - plenum_x_lower

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

formatter = (4*"{:>12}")

print(formatter.format('', 'ncells', 'length', 'delta_x'))
print(formatter.format('Tube_x', tube_n_cells, tube_domain_length, round(tube_dx, 5)))
print(formatter.format('plenum_x', plenum_x_n_cells, plenum_domain_length, round(plenum_dx, 5)))
print(formatter.format('plenum_y', plenum_y_n_cells, plenum_y_length, round(plenum_dy, 5)))
print(formatter.format('plenum_z', plenum_z_n_cells, plenum_z_length, round(plenum_dz, 5)))
print("total cells = {}".format(tube_n_cells+ (plenum_x_n_cells*plenum_y_n_cells)))

RunOptions = {
  'cfl': 0.1125,# / float(tube_n_cells / 64),
  'final_time': 0.2,
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

R = 1.0
gamma = 1.4

Equation = {
  'Rspec': R,
  'gamma': gamma
}

# R_ref = 287.4
# p_ref = 10_000.
# T_ref = 300.
# L_ref = 1.0
# rho_ref = p_ref / T_ref / R_ref
# u_ref = math.sqrt(p_ref / rho_ref)
# t_ref = L_ref / u_ref
# # ud->Msq =  u_ref*u_ref / (R_gas*T_ref);

# p0 = 2.0
# rho0 = math.pow(p0, 1.0 / gamma)
# T0 = p0 / rho0
# p = 0.95 * p0
# # T = T0 + ArrheniusKinetics['Q'] * (gamma - 1.0)
# T = 11.290743302923245
# rho = p / T

checkpoint = ''


def Area(xi):
  A0  = 1.0
  A1  = 4.0 # Reference: 3.0; best: 4.0 
  xi0 = 0.0625 - 1.0
  xi1 = 0.75 - 1.0   # Reference: 0.5; best: 0.75
  Ax = 1.0 if xi < xi0 else A0 + (A1-A0)*(xi-xi0)/(xi1-xi0) if xi < xi1 else A1
  return Ax

ControlOptions = {
  # Tube surface
  'surface_area_tube_inlet': (n_tubes * Area(-1.0) * D) / D,
  'surface_area_tube_outlet': (n_tubes * Area(0.0) * D) / D,
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
  'initial_conditions': {
    'initially_filled_x': 0.4,
    'rho_left': 0.125, # should be equal to rho_right from tube, same for pressure
    'p_left': 0.1,
    'rho_right': 0.125,
    'p_right': 0.1
  },
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
    'r_start': 4.0 * r_tube, # +plenum_dy/100.,
    'r_end': 4.0 * r_tube, # +plenum_dy/100.,
    'y_0': y_0,
    'height': inlet_length,
    'angle':  angle # 60. / 180. * math.pi
  } for y_0 in y0s],
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
  # 'massflow_correlation': 0.06 * 4.0,
}


# get Mirrorbox for each Tube where they are connected with the plenum
def PlenumMirrorBox(y0):
  return BoxWhichContains(DomainAroundPoint(GetCenterPoint(plenum_x_lower, y0), [0.0, -1.0], [-0.9 * plenum_x_lower, +1.0]))


Tube_FluxMethod = FluxMethod
Tube_FluxMethod['area_variation'] = Area

Tubes = [{
  'checkpoint': checkpoint if checkpoint == '' else '{}/Tube_{}'.format(checkpoint, i),
  # 'buffer': 0.06,
  'initial_conditions': {
    'initially_filled_x': 0.4,
    'rho_left': 1.0, # should be equal to rho_right from tube, same for pressure
    'p_left': 10.0,
    'rho_right': 0.125,
    'p_right': 0.1
  },
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
} for (i, y_0) in enumerate(y0s)]


tube_intervals = 0.005
plenum_intervals = 0.02

Output = { 
  'outputs': [
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
    'directory': '{}/Plotfiles'.format(outputPath),
    'intervals': [0.001],
    #'frequencies': [1]
  }, {
    'type': 'CounterOutput',
    'frequencies': [1000]
  },
  {
    'type': 'Checkpoint',
    'intervals': [1.0],
    'directory': '{}/Checkpoint/'.format(outputPath)
  }
  ]
}
