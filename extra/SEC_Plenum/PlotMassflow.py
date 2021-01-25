import os
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn')

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

#plt.style.use('seaborn')
tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "font.size": 11,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 11,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7
}
plt.rcParams.update(tex_fonts)

modes = ['cellwise', 'average_mirror_cells', 'average_ghost_cells', 'average_massflow']
boundaries = ['TurbineMassflowBoundaries', 'TurbineMassflowBoundaries_Jirasek']
for boundary in boundaries:
  f_rhos = []
  f_rho_requireds = []
  times_f_rhos = []
  times_f_rho_requireds = []
  for mode in modes:
    #path_to_log = '/home/zenkechr/FVS_develop/FiniteVolumeSolver/build_2D-Release/SEC_Plenum/0000-{}-{}.log'.format(mode, boundary)
    path_to_log = '/srv/public/Maikel/FiniteVolumeSolver/build_2D-Release/SEC_Plenum/0000-{}-{}.log'.format(mode, boundary)
    output_path = '.'
    path_to_f_rho = '{}/F_rho.log'.format(output_path)
    path_to_hom = '{}/Homogenous.log'.format(output_path)
    os.system('grep "F_rho" {} > {}'.format(path_to_log, path_to_f_rho))
    os.system('grep "Required" {} > {}'.format(path_to_log, path_to_hom))

    with open(path_to_f_rho) as f:
      f_rho = np.array([float(line.split()[-1]) for line in f])
    with open(path_to_f_rho) as f:
      times_f_rho = np.array([float(line.split()[-8][:-2]) for line in f])

    with open(path_to_hom) as f:
      f_rho_required = np.array([float(line.split()[-1]) for line in f][1:])
    with open(path_to_hom) as f:
      times_f_rho_required = np.array([float(line.split()[-4][:-2]) for line in f][1:][::2])

    r_tube = 0.015
    D = 2.0 * r_tube
    scale = 1.0 / (3.0 * D)

    f_rho = scale * f_rho
    f_rho_required = np.reshape(f_rho_required, (f_rho.shape[0], 2))
    f_rho_required = f_rho_required[:,0]
    
    f_rhos.append(f_rho)
    f_rho_requireds.append(f_rho_required)
    times_f_rhos.append(times_f_rho)
    times_f_rho_requireds.append(times_f_rho_required)
    # def dist(x, y):
    #   dx = x - y
    #   return np.sqrt(np.sum(dx * dx)) / dx.shape[0]

    # print(np.max(np.abs(f_rho - f_rho_required)))
    # print(np.sum(np.abs(f_rho - f_rho_required)) / f_rho.shape[0])
    # print(dist(f_rho, f_rho_required))
    # print(np.min(f_rho))

  mode_names = ["cell-wise", "average mirror", "average ghost", "average mass flux"]
  f, axs = plt.subplots(nrows=2, ncols=2, figsize=set_size('thesis', fraction=0.95, subplots=(2,2)), sharex=True, sharey=True)
  for ax, f_rho, f_rho_required, times_f_rho, times_f_rho_required, mode in zip(axs.flat, f_rhos, f_rho_requireds, times_f_rhos, times_f_rho_requireds, mode_names):
    l1 = ax.plot(times_f_rho, f_rho, label='real massflow')
    l2 = ax.plot(times_f_rho_required, f_rho_required, label='required massflow')
    ax.set_xlabel('time')
    ax.set_title(mode)
  for ax in axs.flat:
    ax.label_outer()
  lines = l1+l2
  labs=[ l.get_label() for l in lines ]
  f.legend( lines, labs, loc='lower center', bbox_to_anchor=(0.5,-0.05), 
               bbox_transform=plt.gcf().transFigure, ncol=2, fancybox=True, shadow=True)
  
  # axs[1,0].legend(lines, labs, loc='upper center', 
            #  bbox_to_anchor=(1.1, -0.15),fancybox=False, shadow=False, ncol=3)
  plt.tight_layout()
  plt.savefig('{}/plot-{}.png'.format(output_path, boundary), bbox_inches = 'tight')
  plt.savefig('{}/plot-{}.pdf'.format(output_path, boundary), bbox_inches = 'tight')
  plt.clf()