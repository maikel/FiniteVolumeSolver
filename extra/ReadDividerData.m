function [X, Y, time, data] = ReadDividerData(directory_path)
  files_structure = dir(sprintf('%s%s', directory_path, '/*.dat'));
  ntimesteps = length(files_structure);
  if (ntimesteps == 0) 
     ME = MException('ReadTubeData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  filename = sprintf('%s/%s', directory_path, all_filenames{1});
  [X, Y, ~, ~, ~] = ReadTimestepData2(filename);

  nx = length(X);
  ny = length(Y);

%   columns = timestep.colheaders;
  time = zeros(ntimesteps, 1);
  data.rho = zeros(nx, ny, ntimesteps);
  data.rhou = zeros(nx, ny, ntimesteps);
  data.rhov = zeros(nx, ny, ntimesteps);
  data.rhoE = zeros(nx, ny, ntimesteps);
  data.p = zeros(nx, ny, ntimesteps);
  data.speed_of_sound = zeros(nx, ny, ntimesteps);

  for k = 1 : length(all_filenames)
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    [~, ~, ~, t, indata] = ReadTimestepData2(filename);
    time(k) = t;
    data.rho(:, :, k) = indata(:, :, 1);
    data.rhou(:, :, k) = indata(:, :, 2);
    data.rhov(:, :, k) = indata(:, :, 3);
    data.rhoE(:, :, k) = indata(:, :, 4);
    data.p(:, :, k) = indata(:, :, 5);
    data.speed_of_sound(:, :, k) = indata(:, :, 6);
  end
end

