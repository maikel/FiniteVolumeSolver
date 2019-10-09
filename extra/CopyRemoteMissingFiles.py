from os.path import isfile, join
from subprocess import check_output

# Configure Paths
local_path = '.'
remote_path = '/group/ag_klima/SFB1029_C01/Bericht/1761444/Matlab/Plenum_y0'

remote_server = 'hlrn'
chunk_size = 50

extensions = ['dat.bin', 'dat']

# Get all Remote Files from Server
remote_files = check_output('ssh -t fub-andorra "ls -1 ' + remote_path + ' | tr \'\n\' \', \'"', shell=True)
remote_files = remote_files.decode("utf-8").split(',')

# List of missing files
missing = []

for f in remote_files:
    # split files by name and extension
    a = f.split('.', 1)

    # continue if file does not have an extension
    if len(a) < 2:
        continue

    # if file extension is in list
    if a[1] in extensions:
        if not isfile(join(local_path, f)):
            missing.append(f);

    # continue if not our extension
    else:
        continue

# exit if there are no files missing
if len(missing) == 0:
    print("No files missing!")
    exit()

# ask to copy from remote
decision = input(f"{len(missing)} files are missing. Copy from {remote_server}? (yes/no): ")
if decision == "yes":
    # chunk missing files into 50 each iteration
    for i in range(0, len(missing), chunk_size):
        # join filenames and remote_path and convert to list for scp
        chunk = [join(remote_path, m) for m in missing[i:i+chunk_size]]
        copyfiles = ' '.join(chunk)
        print(f"Copying Chunk {int(i/chunk_size)+1} / {int(len(missing)/chunk_size)+1}")
        check_output('scp ' + remote_server + ':"' + copyfiles + '" "' + local_path + '"', shell=True)
else:
    print("No action performed.")