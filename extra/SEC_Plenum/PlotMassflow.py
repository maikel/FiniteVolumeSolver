import os
import numpy as np
import matplotlib.pyplot as plt

modes = ['cellwise', 'average_mirror_cells', 'average_ghost_cells', 'average_massflow']
boundaries = ['TurbineMassflowBoundaries', 'TurbineMassflowBoundaries_Jirasek']
for boundary in boundaries:
  for mode in modes:
    path_to_log = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Release/SEC_Plenum/0000-{}-{}.log'.format(mode, boundary)
    output_path = '.'
    path_to_f_rho = '{}/F_rho.log'.format(output_path)
    path_to_hom = '{}/Homogenous.log'.format(output_path)

    os.system('grep "F_rho" {} > {}'.format(path_to_log, path_to_f_rho))
    os.system('grep "Required" {} > {}'.format(path_to_log, path_to_hom))

    with open(path_to_f_rho) as f:
      f_rho = np.array([float(line.split()[-1]) for line in f])
      times_f_rho = np.array([float(line.split()[-8][:-2]) for line in f])

    with open(path_to_hom) as f:
      f_rho_required = np.array([float(line.split()[-1]) for line in f][1:])
      times_f_rho_required = np.array([float(line.split()[-4][:-2]) for line in f][1:][::2])

    r_tube = 0.015
    D = 2.0 * r_tube
    scale = 1.0 / (3.0 * D)

    f_rho = scale * f_rho
    f_rho_required = np.reshape(f_rho_required, (f_rho.shape[0], 2))
    f_rho_required = f_rho_required[:,0]
    print(f_rho)
    print(f_rho_required)

    def dist(x, y):
      dx = x - y
      return np.sqrt(np.sum(dx * dx)) / dx.shape[0]

    print(np.max(np.abs(f_rho - f_rho_required)))
    print(np.sum(np.abs(f_rho - f_rho_required)) / f_rho.shape[0])
    print(dist(f_rho, f_rho_required))
    print(np.min(f_rho))

    plt.plot(times_f_rho, f_rho, label='real massflow')
    plt.plot(times_f_rho_required, f_rho_required, label='required massflow')
    plt.xlabel('time')
    plt.ylabel('massflow')
    plt.legend()
    plt.title(mode)
    plt.savefig('{}/plot_{}-{}.png'.format(output_path, mode, boundary))
    plt.clf()