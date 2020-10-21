import math

RunOptions = {
  'cfl': 0.8,
  # 'final_time': 0.2 / math.sqrt(101325.),
  'final_time': 0.2,
}

GridGeometry = {
  'cell_dimensions': [100, 1, 1],
  'coordinates': {
    'lower': [-0.50, -0.015, -0.015],
    'upper': [+0.50, +0.015, +0.015],
  },
  'periodicity': [0, 0, 0]
}

PatchHierarchy = {
 'max_number_of_levels': 1,
 'refine_ratio': [2, 1, 1],
 'blocking_factor': [8, 1, 1],
 'max_grid_size': [200, 1, 1],
}

# reconstruction = "Characteristics"
# reconstruction = "Primitive"
reconstruction = "ConservativeM"
# reconstruction = "Conservative"
# reconstruction = "HLLEM"
# reconstruction = "HLLE"
# reconstruction = "PerfectGas"

paths = {
  'HLLE': './HLLE.h5',
  'HLLEM': './HLLEM.h5',
  'Conservative': './Conservative.h5',
  'ConservativeM': './ConservativeM.h5',
  'Primitive': './IdealGasMix_Prim.h5',
  'Characteristics': './Characteristics3.h5',
  'PerfectGas': './PerfectGas1d_cfl05.h5'
}

Output = {
  'outputs': [{
    'type': 'HDF5',
    'path': paths[reconstruction],
    'intervals': [RunOptions['final_time'] / 10.0],
    # 'frequencies': [1],
  }]
}
