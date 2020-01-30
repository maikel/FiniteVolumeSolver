BK19LevelIntegrator = {
    'mlmg_tolerance_rel': 1e-5,
    'mlmg_tolerance_abs': 1e-5,
    'mlmg_verbose': 1,
    'mlmg_max_iter': 1,
    'bottom_tolerance_rel': 1e-7,
    'bottom_tolerance_abs': -1.0,
    'bottom_verbose': 1,
    'bottom_max_iter': 300
}

GridGeometry = {
    'cell_dimensions': [32, 32],
#    'cell_dimensions': [64, 64],
    'x_lower': [0.0, 0.0],
    'x_upper': [1.0, 1.0],
    'periodicity': [1, 1]
}

PatchHierarchyOptions = {
    'max_number_of_levels': 1,
    'refine_ratio': [2, 2],
    'blocking_factor': [32, 32],
#    'max_grid_size': [64, 64]
    'max_grid_size': [32, 32]
}
