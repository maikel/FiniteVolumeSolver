function [X, Y, time, data, columns] = ReadPlenumData(directory_path)
  files_structure = dir(sprintf('%s%s', directory_path, '/*.dat'));
  ntimesteps = length(files_structure);
  if (ntimesteps == 0) 
     ME = MException('ReadPlenumData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  timestep = importdata(sprintf('%s/%s', directory_path, all_filenames{1}), ' ', 5);
  nx_string = timestep.textdata{1};
  nx = sscanf(nx_string, 'nx = %d');
  ny_string = timestep.textdata{2};
  ny = sscanf(ny_string, 'ny = %d');
  nspecies = length(timestep.colheaders) - 11;
  
  columns = timestep.colheaders;
  
  X = reshape(timestep.data(:, 1), nx, ny);
  Y = reshape(timestep.data(:, 2), nx, ny);
  time = zeros(ntimesteps, 1);
  data.rho = zeros(nx, ny, ntimesteps);
  data.u = zeros(nx, ny, ntimesteps, 3);
  data.p = zeros(nx, ny, ntimesteps);
  data.speed_of_sound = zeros(nx, ny, ntimesteps);
  data.T = zeros(nx, ny, ntimesteps);
  data.gamma = zeros(nx, ny, ntimesteps);
  data.heat_capacity_at_constant_pressure = zeros(nx, ny, ntimesteps);
  data.species = zeros(nx, ny, ntimesteps, nspecies);

  for k = 1 : ntimesteps
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    timestep = importdata(filename, ' ', 5);
    timepoint_string = timestep.textdata{3};
    time(k) = sscanf(timepoint_string, 't = %f');
    data.rho(:, :, k) = reshape(timestep.data(:, 3), nx, ny);
    data.u(:, :, k, 1) = reshape(timestep.data(:, 4), nx, ny);
    data.u(:, :, k, 2) = reshape(timestep.data(:, 5), nx, ny);
    data.u(:, :, k, 3) = reshape(timestep.data(:, 6), nx, ny);
    data.speed_of_sound(:, :, k) = reshape(timestep.data(:, 7), nx, ny);
    data.T(:, :, k) = reshape(timestep.data(:, 8), nx, ny);
    data.p(:, :, k) = reshape(timestep.data(:, 9), nx, ny);
    data.gamma(:, :, k) = reshape(timestep.data(:, 10), nx, ny);
    data.heat_capacity_at_constant_pressure(:, :, k) = reshape(timestep.data(:, 11), nx, ny);
    data.species(:, :, k, :) = reshape(timestep.data(:, 12:end), nx, ny, nspecies);
  end
end
