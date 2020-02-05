import math
import h5py
import numpy as np

source = h5py.File("/Users/maikel/Development/FiniteVolumeSolver/build/initial_data2.h5", "r")
dataset = source['DG_Solution']
data = np.transpose(np.array(dataset)).copy()
source.close()

view = data.view()
view.shape = (data.shape[0] * data.shape[1])
data = view.copy()

plenum_y_n_cells = 64
plenum_x_n_cells = 64
plenum_z_n_cells = 64

tube_blocking_factor = 8
plenum_blocking_factor = 8

n_level = 1

n_tubes = 6
r_tube = 0.015
r_inner = 0.5 * 0.130
r_outer = 0.5 * 0.385
r_tube_center = 0.5 * r_inner + 0.5 * r_outer
alpha = 2.0 * math.pi / n_tubes

tube_length = 1.50 # [m]
plenum_length = 0.50 # [m]

plenum_max_grid_size = max(plenum_blocking_factor, 64)

plenum_domain_length = plenum_length + tube_length
plenum_x_upper = 0.25 # plenum_length
plenum_x_lower = -0.1 # -tube_length

plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_y_upper = +r_outer + 0.07
plenum_y_lower = 0.0 # -r_outer - 0.05
plenum_y_length = plenum_y_upper - plenum_y_lower

plenum_z_upper = r_outer + 0.05
plenum_z_lower = 0.0
plenum_z_length = plenum_z_upper - plenum_z_lower


plenum_x_over_z_ratio = plenum_x_length / plenum_z_length
plenum_y_over_z_ratio = plenum_y_length / plenum_z_length

plenum_x_n_cells = plenum_z_n_cells * plenum_x_over_z_ratio
plenum_x_n_cells -= plenum_x_n_cells % plenum_blocking_factor
plenum_x_n_cells = int(plenum_x_n_cells)

plenum_y_n_cells = plenum_z_n_cells * plenum_y_over_z_ratio
plenum_y_n_cells -= plenum_y_n_cells % plenum_blocking_factor
plenum_y_n_cells = int(plenum_y_n_cells)


RunOptions = {
  'cfl': 0.3,
  'final_time': 0.02,
  'max_cycles': 10,
  'do_backup': 0
}

# checkpoint = '/Users/maikel/Development/FiniteVolumeSolver/build_3d/MultiTube/Checkpoint/000000010'
checkpoint = ''

GridGeometry = {
  'cell_dimensions': [plenum_x_n_cells, plenum_y_n_cells, plenum_z_n_cells],
  'coordinates': {
    'lower': [plenum_x_lower, plenum_y_lower, plenum_z_lower],
    'upper': [plenum_x_upper, plenum_y_upper, plenum_z_upper],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
  'max_number_of_levels': n_level, 
  'blocking_factor': [plenum_blocking_factor, plenum_blocking_factor, plenum_blocking_factor],
  'max_grid_size': [plenum_max_grid_size, plenum_max_grid_size, plenum_max_grid_size],
  'ngrow_eb_level_set': 5
}

InitialCondition = {
  'data': data
}

PressureBoundary = {
  'outer_pressure': 101325.0,
  'side': 1,
  'direction': 0
}

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'MultiTube/',
    # 'intervals': [1e-4],
    'frequencies': [1],
  }, {
    'type': 'CounterOutput',
    'frequencies': [20]
  }, {
    'type': 'Checkpoint',
    'directory': 'MultiTube/Checkpoint/',
    'frequencies': [100]
  }]
}

