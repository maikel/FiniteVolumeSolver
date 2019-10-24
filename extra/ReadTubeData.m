function [X, time, data] = ReadTubeData(directory_path)
  files_structure = dir(sprintf('%s%s', directory_path, '/*.dat'));
  ntimesteps = length(files_structure);
  if (ntimesteps == 0) 
     ME = MException('ReadTubeData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  timestep = importdata(sprintf('%s/%s', directory_path, all_filenames{1}), ' ', 4);
  nx = size(timestep.data, 1);
  nspecies = size(timestep.data, 2) - 8;

%   columns = timestep.colheaders;
  X = timestep.data(:, 1);
  time = zeros(length(all_filenames), 1);
  data.rho = zeros(nx, length(all_filenames));
  data.u = zeros(nx, length(all_filenames));
  data.p = zeros(nx, length(all_filenames));
  data.T = zeros(nx, length(all_filenames));
  data.speed_of_sound = zeros(nx, ntimesteps);
  data.T = zeros(nx, ntimesteps);
  data.gamma = zeros(nx, ntimesteps);
  data.heat_capacity_at_constant_pressure = zeros(nx, ntimesteps);
  data.species = zeros(nx, ntimesteps, nspecies);

  for k = 1 : length(all_filenames)
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    fprintf('Process %s\n', filename);
    timestep = importdata(filename, ' ', 4);
    timepoint_string = timestep.textdata{2};
    time(k) = sscanf(timepoint_string, 't = %f');
    data.rho(:, k) = timestep.data(:, 2);
    data.u(:, k) = timestep.data(:, 3);
    data.speed_of_sound(:, k) = timestep.data(:, 4);
    data.T(:, k) = timestep.data(:, 5);
    data.p(:, k) = timestep.data(:, 6);
    data.gamma(:, k) = timestep.data(:, 7);
    data.heat_capacity_at_constant_pressure(:, k) = timestep.data(:, 8);
    data.species(:, k, :) = reshape(timestep.data(:, 9:end), nx, nspecies);
  end
end

