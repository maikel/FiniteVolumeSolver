import math
import h5py
import numpy as np

source = "/home/mi/guttula/initial_data2.h5"
#source = h5py.File("/Users/maikel/Development/FiniteVolumeSolver/initial_data2.h5", "r")
#dataset = source['DG_Solution']
#data = np.array([-1.0, 1.0, 1.0, 0.0, 0.0, 1.0])
#data2 = np.array(dataset)
#source.close()

#view = data.view()
#view.shape = (data.shape[0] * data.shape[1])
#data = view.copy()

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
plenum_x_upper = plenum_length
plenum_x_lower = -tube_length

plenum_x_length = plenum_x_upper - plenum_x_lower

plenum_y_upper = +r_outer + 0.005
plenum_y_lower = 0.0 # -r_outer - 0.05
plenum_y_length = plenum_y_upper - plenum_y_lower

plenum_z_upper = r_outer + 0.005
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
  'cfl': 0.5,
  'final_time': 0.02,
  'max_cycles': 1,
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
  'ngrow_eb_level_set': 5,
  'cutcell_load_balance_weight': 10
}

InitialCondition = {
  'data': source
}

PressureBoundary = {
  'outer_pressure': 101325.0,
  'side': 1,
  'direction': 0
}

probe_xs = [
        [ 0.43      ,  0.16789502,  0.09416822],
        [ 0.43      ,  0.17748883,  0.07452529],
        [ 0.43      ,  0.18479092,  0.05393096],
        [ 0.43      ,  0.18971224,  0.03264152],
        [ 0.43      ,  0.19218914,  0.0109269 ],
        [ 0.43      ,  0.19218919, -0.0109261 ],
        [ 0.43      ,  0.18971283, -0.03263811],
        [ 0.43      ,  0.18479114, -0.05393019],
        [ 0.43      ,  0.17748915, -0.07452455],
        [ 0.43      ,  0.16789635, -0.09416584],
        [ 0.43      ,  0.15579194,  0.08737832],
        [ 0.43      ,  0.16469269,  0.06915219],
        [ 0.43      ,  0.17146924,  0.05004325],
        [ 0.43      ,  0.17603584,  0.03028613],
        [ 0.43      ,  0.17833387,  0.01013892],
        [ 0.43      ,  0.17833391, -0.01013817],
        [ 0.43      ,  0.17603597, -0.0302854 ],
        [ 0.43      ,  0.17146946, -0.05004254],
        [ 0.43      ,  0.16469299, -0.06915151],
        [ 0.43      ,  0.15579231, -0.08737765],
        [ 0.43      ,  0.14368758,  0.08058987],
        [ 0.43      ,  0.15189775,  0.0637802 ],
        [ 0.43      ,  0.15814695,  0.04615534],
        [ 0.43      ,  0.16235771,  0.02793383],
        [ 0.43      ,  0.16447881,  0.00935153],
        [ 0.43      ,  0.16447881, -0.00935085],
        [ 0.43      ,  0.16235792, -0.02793311],
        [ 0.43      ,  0.15814716, -0.04615466],
        [ 0.43      ,  0.15189802, -0.06377956],
        [ 0.43      ,  0.14368793, -0.08058925],
        [ 0.43      ,  0.13158376,  0.07380093],
        [ 0.43      ,  0.13910217,  0.05840753],
        [ 0.43      ,  0.14482505,  0.04226751],
        [ 0.43      ,  0.14868127,  0.02558037],
        [ 0.43      ,  0.15062459,  0.00856135],
        [ 0.43      ,  0.15062459, -0.00856077],
        [ 0.43      ,  0.14868128, -0.02557976],
        [ 0.43      ,  0.14482524, -0.04226691],
        [ 0.43      ,  0.13910241, -0.05840695],
        [ 0.43      ,  0.13158409, -0.07380034],
        [ 0.43      ,  0.11947974,  0.06701228],
        [ 0.43      ,  0.12630656,  0.0530349 ],
        [ 0.43      ,  0.13150285,  0.03837961],
        [ 0.43      ,  0.13500521,  0.02322745],
        [ 0.43      ,  0.13677001,  0.00777293],
        [ 0.43      ,  0.13676836, -0.00777534],
        [ 0.43      ,  0.1350053 , -0.02322691],
        [ 0.43      ,  0.13150301, -0.03837906],
        [ 0.43      ,  0.12630679, -0.05303437],
        [ 0.43      ,  0.11948004, -0.06701174],
        [ 0.43      ,  0.10737585,  0.06022344],
        [ 0.43      ,  0.11351093,  0.04766231],
        [ 0.43      ,  0.11818083,  0.0344917 ],
        [ 0.43      ,  0.12132861,  0.02087477],
        [ 0.43      ,  0.12291206,  0.00698791],
        [ 0.43      ,  0.12291338, -0.0069881 ],
        [ 0.43      ,  0.12132871, -0.02087428],
        [ 0.43      ,  0.11818097, -0.03449121],
        [ 0.43      ,  0.11351112, -0.04766184],
        [ 0.43      ,  0.10737613, -0.06022294],
        [ 0.43      ,  0.09527165,  0.05343446],
        [ 0.43      ,  0.10071575,  0.04228936],
        [ 0.43      ,  0.10485872,  0.0306039 ],
        [ 0.43      ,  0.10764911,  0.01852595],
        [ 0.43      ,  0.1090578 ,  0.00619982],
        [ 0.43      ,  0.10905778, -0.0061994 ],
        [ 0.43      ,  0.10764928, -0.01852535],
        [ 0.43      ,  0.10485886, -0.03060347],
        [ 0.43      ,  0.10071591, -0.04228895],
        [ 0.43      ,  0.09527192, -0.05343401],
        [ 0.43      ,  0.08316793,  0.04664586],
        [ 0.43      ,  0.08792025,  0.03691642],
        [ 0.43      ,  0.09153667,  0.0267161 ],
        [ 0.43      ,  0.09397578,  0.01616966],
        [ 0.43      ,  0.09520141,  0.00541324],
        [ 0.43      ,  0.09520143, -0.00541285],
        [ 0.43      ,  0.09397585, -0.01616929],
        [ 0.43      ,  0.09153678, -0.02671571],
        [ 0.43      ,  0.08792046, -0.03691598],
        [ 0.43      ,  0.08316815, -0.04664543],
        [ 0.43      ,  0.07106422,  0.03985637],
        [ 0.43      ,  0.07512392,  0.0315449 ],
        [ 0.43      ,  0.07821486,  0.02282857],
        [ 0.43      ,  0.08029642,  0.01381815],
        [ 0.43      ,  0.08134657,  0.0046243 ],
        [ 0.43      ,  0.08134656, -0.00462403],
        [ 0.43      ,  0.08029654, -0.01381774],
        [ 0.43      ,  0.0782149 , -0.02282815],
        [ 0.43      ,  0.07512398, -0.0315446 ],
        [ 0.43      ,  0.07106443, -0.039856  ],
        [ 0.43      ,  0.05896053,  0.03306888],
        [ 0.43      ,  0.06232846,  0.02617178],
        [ 0.43      ,  0.0648925 ,  0.01894035],
        [ 0.43      ,  0.06662146,  0.0114621 ],
        [ 0.43      ,  0.06749001,  0.00383642],
        [ 0.43      ,  0.06748999, -0.00383603],
        [ 0.43      ,  0.06662147, -0.01146183],
        [ 0.43      ,  0.06489255, -0.01894009],
        [ 0.43      ,  0.06232857, -0.02617153],
        [ 0.43      ,  0.05896072, -0.03306856]]

probes = filter(lambda x: plenum_z_lower < x[2] and x[2] < plenum_z_upper, probe_xs)
probe_tuples = [(x[0], x[1], x[2]) for x in probes]

Output = { 
  'outputs': [{
    'type': 'Plotfile',
    'directory': 'MultiTube/Plotfile/',
    'intervals': [5e-5]
  }, {
    'type': 'CounterOutput',
    'frequencies': [100]
  },
  {
    'type': 'Checkpoint',
    'directory': 'MultiTube/Checkpoint/',
    'intervals': [3e-4]
  }, {
     'type': 'ProbesOutput',
     'filename': 'Probes.h5',
     'probes': probe_tuples,
     'frequencies': [2]
  }]
}

