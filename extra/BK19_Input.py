BK19LevelIntegrator = {
    'mlmg_tolerance_rel': 0.001,
    'mlmg_tolerance_abs': -1.0,
    'mlmg_verbose': 4,
    'mlmg_max_iter': 1,
    'bottom_tolerance_rel': 1e-4,
    'bottom_tolerance_abs': -1.0,
    'bottom_verbose': 4,
    'bottom_max_iter': 1000
}

GridGeometry = {
    'cell_dimensions': [64, 64],
    'x_lower': [0.0, 0.0],
    'x_upper': [1.0, 1.0],
    'periodicity': [1, 1]
}

PatchHierarchyOptions = {
    'max_number_of_levels': 1,
    'refine_ratio': [2, 2],
    'blocking_factor': [64, 64],
    'max_grid_size': [64, 64]
}
