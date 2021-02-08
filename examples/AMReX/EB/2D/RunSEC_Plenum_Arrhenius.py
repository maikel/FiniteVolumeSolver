# /usr/bin/python

import os, sys
import shutil

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'

os.environ['OMP_NUM_THREADS'] = '1'

src_dir = FVS_path
build_dir = "{}/build_2D-Release".format(src_dir)
work_dir = "{}/SEC_Plenum_Arrhenius".format(build_dir)

application_dir = "{}/examples".format(build_dir)
application_src = "{}/AMReX.EB.SEC_Plenum_Arrhenius".format(application_dir)
application_dest = "{}/AMReX.EB.SEC_Plenum_Arrhenius".format(work_dir)

input_file = "{}/examples/AMReX/EB/2D/SEC_Plenum_Arrhenius.py".format(src_dir)

os.makedirs(work_dir, exist_ok=True)
os.chdir(work_dir)
shutil.copy(application_src, work_dir)

mpi_command = 'mpiexec -n 6'

mul = [3.0, 1.0, 0.5, 0.1]
# mul = [0.1]

massflow = [0.4, 0.6, 1.0]
# massflow = [0.1]

for mu in mul:
  for mass in massflow:
    relative_outputpath = 'mul_{}_massflow_{}'.format(str(mu).replace('.', '-'), str(mass).replace('.', '-') )
    os.makedirs(relative_outputpath, exist_ok=True)

    input_file_workdir = '{}/{}/SEC_Plenum_Arrhenius.py'.format(work_dir, relative_outputpath )
    with open(input_file) as f:
      newText = f.read().replace('%Diffusion%', str(mu))
      newText = newText.replace('%MASSFLOW%', str(mass))
      newText = newText.replace('%OUTPUT%', str(relative_outputpath))
    
    with open(input_file_workdir, "w") as f:
      f.write(newText)
    os.system('{} {} --config {}'.format(mpi_command, application_dest, input_file_workdir))
    os.system('cp 0000.log {}/0000-mu_{}_mass_{}.log'.format(relative_outputpath, str(mu).replace('.', '-'), str(mass).replace('.', '-') ) )
    os.system('mv {} {}'.format('Checkpoint', relative_outputpath) )

for mu in mul:
  for mass in massflow:
    relative_outputpath = 'mul_{}_massflow_{}'.format(str(mu).replace('.', '-'), str(mass).replace('.', '-') )

    os.chdir(FVS_path+"/extra/SEC_Plenum/")

    os.system('python3 PlotControlState_HDF5.py {}/{}/'.format(work_dir, relative_outputpath) )
    os.system('python3 PlotTube_HDF5.py {}/{}/'.format(work_dir, relative_outputpath) )
    # os.system('python3 PlotPlenum_HDF5.py {}/{}/'.format(work_dir, relative_outputpath) )