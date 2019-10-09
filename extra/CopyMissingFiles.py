from os import listdir
from os.path import isfile, join
from subprocess import check_output

# Configure Paths
local_path='.'
remote_path='/scratch/usr/beinadol/Turbine_WithOffset/1761444/MultiTube/Matlab/Plenum_y0/'

remote_server = 'hlrn'
chunk_size = 50

# Get local files
local_files = [f for f in listdir(local_path) if isfile(join(local_path, f))]

# List of missing files
missing = []

for f in local_files:
    # split files by name and extension
    a = f.split('.', 1)

    # continue if file does not have an extension
    if len(a) < 2:
        continue

    # if file extension is .dat.bin then look for corresponding .dat file
    if a[1] == 'dat.bin':
        if not isfile(join(local_path, a[0]+'.dat')):
            missing.append(a[0]+'.dat');

    # if file extension is .dat then look for corresponding .dat.bin file
    elif a[1] == 'dat':
        if not isfile(join(local_path, a[0]+'.dat.bin')):
            missing.append(a[0]+'.dat.bin')

    # continue if not our extension
    else:
        continue

# exit if there are no files missing
if len(missing) == 0:
    print("No files missing!")
    exit()

# ask to copy from remote
decision = input("{} files are missing. Copy from {}? (yes/no): ".format(len(missing), remote_server))
if decision == "yes":
    # chunk missing files into 50 each iteration
    for i in range(0, len(missing), chunk_size):
        # join filenames and remote_path and convert to list for scp
        chunk = [join(remote_path, m) for m in missing[i:i+chunk_size]]
        copyfiles = ' '.join(chunk)
        print("Copying Chunk {} / {}".format(int(i/chunk_size)+1, int(len(missing)/chunk_size)+1))
        check_output('scp ' + remote_server + ':"' + copyfiles + '" "' + local_path + '"', shell=True)
else:
    print("No action performed.")
