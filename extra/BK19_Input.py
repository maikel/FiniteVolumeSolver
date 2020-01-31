BK19LevelIntegrator = {
    'mlmg_tolerance_rel': -1,
    'mlmg_tolerance_abs': 1e-12,
    'mlmg_verbose': 0,
    'mlmg_max_iter': 100,
    'bottom_tolerance_rel': 1e-6,
    'bottom_tolerance_abs': -1.0,
    'bottom_verbose': 0,
    'bottom_max_iter': 50
}

GridGeometry = {
    # 'cell_dimensions': [32, 32],
    'cell_dimensions': [32, 32],
    'x_lower': [0.0, 0.0],
    'x_upper': [1.0, 1.0],
    'periodicity': [1, 1]
}

PatchHierarchy = {
    'max_number_of_levels': 1,
    'refine_ratio': [2, 2],
    'blocking_factor': [32, 32],
    'max_grid_size': [256, 256]
}

RunOptions = {
    'cfl': 0.45,
    'final_time': 1.0,
}

Output = {
    'outputs': [
        {
            'type': 'Plotfile',
            'intervals': [5e-2]
        },
    ]
}