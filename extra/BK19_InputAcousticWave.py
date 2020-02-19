BK19LevelIntegrator = {
    'mlmg_tolerance_rel': 1e-10,
    'mlmg_tolerance_abs': 1e-10,
    'mlmg_verbose': 2,
    'mlmg_max_iter': 100,
    'bottom_tolerance_rel': 1e-4,
    'bottom_tolerance_abs': -1.0,
    'bottom_verbose': 2,
    'bottom_max_iter': 200
}

GridGeometry = {
    'cell_dimensions': [64, 64],
    #'cell_dimensions': [32, 32],
    'x_lower': [0.0, 0.0],
    'x_upper': [1.0, 1.0],
    'periodicity': [1, 1]
}

PatchHierarchy = {
    'max_number_of_levels': 1,
    'refine_ratio': [2, 2],
    #'blocking_factor': [64, 64],
    'blocking_factor': [32, 32],
    'max_grid_size': [256, 256]
}

RunOptions = {
    'cfl': 0.45,
    'final_time': 0.002,
}

Output = {
    'outputs': [
        {
            'type': 'Plotfile',
            #'intervals': [5e-5]
            'intervals': [2.5e-05]
        },
    ]
}
