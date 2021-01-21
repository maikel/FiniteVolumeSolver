# /usr/bin/python

import os
import shutil

# export OMP_NUM_THREADS=1
os.environ['OMP_NUM_THREADS'] = '1'

src_dir = "/srv/public/Maikel/FiniteVolumeSolver"
build_dir = "{}/build_2D-Release".format(src_dir)
work_dir = "{}/SEC_Plenum".format(build_dir)

application_dir = "{}/examples".format(build_dir)
application_src = "{}/AMReX.EB.SEC_Plenum".format(application_dir)
application_dest = "{}/AMReX.EB.SEC_Plenum".format(work_dir)

input_file = "{}/examples/AMReX/EB/2D/SEC_Plenum.py".format(src_dir)

os.makedirs(work_dir, exist_ok=True)
os.chdir(work_dir)
shutil.copy(application_src, work_dir)

mpi_command = 'mpiexec'
modes = [0, 1, 2, 3]
mode_names = ['cellwise', 'average_mirror_cells', 'average_ghost_cells', 'average_massflow']
boundaries = ['TurbineMassflowBoundaries', 'TurbineMassflowBoundaries_Jirasek']

for boundary in boundaries:
  for mode in modes:
      input_file_workdir = '{}/SEC_Plenum_{}-.py'.format(work_dir, mode_names[mode], boundary)
      with open(input_file) as f:
        newText = f.read().replace('%MODE%', str(mode)).replace('%BOUNDARY_CONDITION%', boundary)
      with open(input_file_workdir, "w") as f:
        f.write(newText)
      os.system('{} {} --config {}'.format(mpi_command, application_dest, input_file_workdir))
      os.system('cp 0000.log 0000-{}-{}.log'.format(mode_names[mode], boundary))
