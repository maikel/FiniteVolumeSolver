function [X, Y, time, data] = ReadDividerData(directory_path)
  files_structure = dir(sprintf('%s%s', directory_path, '/*.dat'));
  ntimesteps = length(files_structure);
  if (ntimesteps == 0) 
     ME = MException('ReadTubeData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  timestep = importdata(sprintf('%s/%s', directory_path, all_filenames{1}), ' ', 5);
  nx = sscanf(timestep.textdata{1},'nx = %d');
  ny = sscanf(timestep.textdata{2,1},'ny = %d');

%   columns = timestep.colheaders;
  indata = reshape(timestep.data, nx, ny, 7);
  X = indata(:, 1, 1);
  Y = indata(1, :, 2);
  time = zeros(length(all_filenames), 1);
  data.rho = zeros(nx, ny, length(all_filenames));
  data.u = zeros(nx, ny, length(all_filenames));
  data.v = zeros(nx, ny, length(all_filenames));
  data.p = zeros(nx, ny, length(all_filenames));
  data.speed_of_sound = zeros(nx, ny, ntimesteps);

  for k = 1 : length(all_filenames)
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    timestep = importdata(filename, ' ', 5);
    indata = reshape(timestep.data, nx, ny, 7);
    timepoint_string = timestep.textdata{3};
    time(k) = sscanf(timepoint_string, 't = %f');
    data.rho(:, :, k) = indata(:, :, 3);
    data.u(:, :, k) = indata(:, :, 4);
    data.v(:, :, k) = indata(:, :, 5);
    data.speed_of_sound(:, :, k) = indata(:, :, 6);
    data.p(:, :, k) = indata(:, :, 7);
  end
end

