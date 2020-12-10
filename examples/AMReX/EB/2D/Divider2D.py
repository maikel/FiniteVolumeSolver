import math

def nCellsY(nCellsX, ratio):
    base = int(nCellsX * ratio)
    ncells = base - base % 8
    return ncells

cfl = 0.4
final_time = 1.01e-3

Mach_number = 1.1

xlower = 0.0
xupper = 0.2300001
xlen = xupper - xlower

ylower = -0.05 + 0.005
yupper = +0.05 + 0.005
ylen = yupper - ylower

nx = 800
ny = nCellsY(nx, ylen / xlen)

n_cells = [nx, ny]
n_level = 1
x_range = [xlower, xupper]
y_range = [ylower, yupper]

base_path = '/srv/public/Maikel/FiniteVolumeSolver/examples/AMReX/EB/2D/'

wall_filenames = ['{}/wall_1.txt'.format(base_path),
                  '{}/wall_2.txt'.format(base_path),
                  '{}/wall_3.txt'.format(base_path)]

# max_cycles = 1

output = {
	'outputs': [
  # {'type': 'Checkpoint', 'directory': 'Divider_DE5/Checkpoint', 'intervals': [1e-4]},
  # {'type': 'Plotfiles', 'directory': 'Divider_DE5/Plotfiles', 'intervals': [1e-4]},
  {'type': 'HDF5', 'path': 'Divider_c24.h5', 'intervals': [1e-5]}
]
}
