BK19Solver = {
    'do_initial_projection': True,
    'mlmg_tolerance_rel': 1e-10,
    'mlmg_tolerance_abs': 1e-10,
    'mlmg_verbose': 0,
    'mlmg_max_iter': 100,
    'bottom_tolerance_rel': 1e-4,
    'bottom_tolerance_abs': -1.0,
    'bottom_verbose': 0,
    'bottom_max_iter': 200
}

R_gas = 287.0
h_ref = 1.0
t_ref = 1.0
T_ref = 353.048780488
u_ref = h_ref / t_ref

gamma = 2.0

BK19PhysicalParameters = {
    'R_gas': R_gas,
    'gamma': gamma,
    'Msq': u_ref * u_ref / R_gas * T_ref,
    'c_p': gamma / (gamma - 1.0),
    'alpha_p': 1.0
}

LinearOperator = 'MLNodeHelmDualCstVel'

GridGeometry = {
    'cell_dimensions': [8, 8, 8],
    #'cell_dimensions': [64, 64, 64],
    #'cell_dimensions': [32, 32, 32],
    'coordinates': {
      'lower': [0.0, 0.0, 0.0],
      'upper': [1.0, 1.0, 1.0]
    },
    'periodicity': [1, 0, 1]
}

PatchHierarchy = {
    'max_number_of_levels': 1,
    'refine_ratio': [2, 2, 2],
    #'blocking_factor': [64, 64, 64],
    'blocking_factor': [32, 32, 32],
    'max_grid_size': [256, 256, 256]
}

RunOptions = {
    'cfl': 0.45,
    'final_time': 0.002
}

Output = {
    'outputs': [
        {
            'type': 'DebugOutput',
            'directory': './AcousticWaveWall/Debug/'
        },
        {
            'type': 'Plotfile',
            'directory': './AcousticWaveWall/Plotfile/',
            #'intervals': [5e-5]
            'intervals': [2.5e-05]
        },
    ]
}
