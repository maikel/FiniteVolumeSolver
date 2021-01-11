import math

plenum_x_n_cells = 128
tube_blocking_factor = 8
plenum_blocking_factor = 8

n_level = 1

n_tubes = 6
r_tube = 0.015

D = 2.0 * r_tube

r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.2
r_tube_center = 2.0 * r_inner
alpha = 2.0 * math.pi / n_tubes

inlet_length = 3.0 * D # [m]

plenum_x_upper = 0.1
plenum_x_lower = -inlet_length
plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_length = plenum_x_length # [m]
tube_length = 1.0 # [m]

plenum_max_grid_size = max(plenum_blocking_factor, 1024)

plenum_domain_length = plenum_length + inlet_length
tube_domain_length = tube_length - inlet_length

tube_over_plenum_length_ratio = tube_domain_length / plenum_domain_length

# plenum_yz_upper = +(r_outer + 0.01)
# plenum_yz_lower = -plenum_yz_upper

# plenum_yz_length = plenum_yz_upper - plenum_yz_lower


plenum_y_lower = - 0.5
plenum_y_upper = + 0.5
plenum_y_length = plenum_y_upper - plenum_y_lower

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

tube_n_cells = plenum_x_n_cells * tube_over_plenum_length_ratio
tube_n_cells -= tube_n_cells % tube_blocking_factor
tube_n_cells = int(tube_n_cells)

RunOptions = {
  'cfl': 0.8,
  'final_time': 20,
  'max_cycles': -1
}

# checkpoint = '/Users/maikel/Development/FiniteVolumeSolver/build_3d/MultiTube/Checkpoint/000000063'
checkpoint = ''

def ToCellIndex(x, xlo, xhi, ncells):
  xlen = xhi - xlo
  x_rel = (x - xlo) / xlen
  i = int(x_rel * ncells)
  return i

y0s = [-1.0/3.0, 0.0, +1.0/3.0]
mach_1_boundaries = [ToCellIndex(y0, plenum_y_lower, plenum_y_upper, plenum_y_n_cells) for y0 in y0s]

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
  'FluxMethod': {
    'limiter': 'Upwind',
  },
  'InletGeometries': [{
    'r_start': r_tube,
    'r_end': 2.0*r_tube,
    'y_0': y_0,
  } for y_0 in y0s],
  'InitialCondition': {
    'left': {
      'density': 1.0,
      'temperature': 1.0
    },
    'right': {
      'density': 1.0,
      'temperature': 1.0
    },
  },
  #'PressureOutflowBoundaries': [{
  #  'outer_pressure': 1.0,
  #  'efficiency': 0.1
  #}],
  'TurbineMassflowBoundaries': [{
    'boundary_section': { 
      'lower': [plenum_x_n_cells, -4, 0], 
      'upper': [plenum_x_n_cells + 3, plenum_y_n_cells + 3, 0] 
     },
    'relative_surface_area': 1,
    'massflow_correlation': 0.02,
  }],
  #'MachnumberBoundaries': [{
  #  'boundary_section': { 
  #    'lower': [plenum_x_n_cells, y0 - int(r_tube / plenum_y_length * plenum_y_n_cells), 0], 
  #    'upper': [plenum_x_n_cells + 1, y0 + int(r_tube / plenum_y_length * plenum_y_n_cells), 0] 
  #   }
  #} for y0 in mach_1_boundaries]
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
  print(real_box)
  i0 = ToCellIndex(real_box[0][0], plenum_x_lower, plenum_x_upper, plenum_x_n_cells)
  iEnd = ToCellIndex(real_box[1][0], plenum_x_lower, plenum_x_upper, plenum_x_n_cells)
  j0 = ToCellIndex(real_box[0][1], plenum_y_lower, plenum_y_upper, plenum_y_n_cells)
  jEnd = ToCellIndex(real_box[1][1], plenum_y_lower, plenum_y_upper, plenum_y_n_cells)
  return { 'lower': [i0, j0, 0], 'upper': [iEnd, jEnd, 0] }

def PlenumMirrorBox(y0):
  return BoxWhichContains(DomainAroundPoint(TubeCenterPoint(-inlet_length, y0), [0.0, -D], [inlet_length, D]))

Tubes = [{
  'checkpoint': checkpoint,
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
  }
} for y_0 in y0s]

Output = { 
  'outputs': [{
    'type': 'Plotfiles',
    'directory': 'SEC_Plenum_TurbineBoundary/Plotfiles/',
    'intervals': [0.01],
    #'frequencies': [1]
  }, {
    'type': 'Checkpoint',
    'directory': 'SEC_Plenum_TurbineBoundary/Checkpoint/',
    'intervals': [1.0],
    #'frequencies': [100]
  }, {
   'type': 'CounterOutput',
   'intervals': [1.0]
  }
  ]
}
