function [X, Y, Z, time, data, columns] = ReadPlenumData(directory_path)
  files_structure = dir(sprintf('%s%s', directory_path, '/*.dat'));
  ntimesteps = length(files_structure);
  if (ntimesteps == 0) 
     ME = MException('ReadPlenumData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  timestep = importdata(sprintf('%s/%s', directory_path, all_filenames{1}), ' ', 6);
  timestep.colheaders = regexp(timestep.textdata{6}, ' +', 'split');
  nx_string = timestep.textdata{1};
  nx = sscanf(nx_string, 'nx = %d');
  ny_string = timestep.textdata{2};
  ny = sscanf(ny_string, 'ny = %d');
  nz_string = timestep.textdata{3};
  nz = sscanf(nz_string, 'nz = %d');
  nspecies = length(timestep.colheaders) - 12;
  
  columns = timestep.colheaders;
  
  X = reshape(timestep.data(:, 1), nx, ny, nz);
  Y = reshape(timestep.data(:, 2), nx, ny, nz);
  Z = reshape(timestep.data(:, 3), nx, ny, nz);
  time = zeros(ntimesteps, 1);
  data.rho = zeros(nx, ny, nz, ntimesteps);
  data.u = zeros(nx, ny, nz, ntimesteps, 3);
  data.p = zeros(nx, ny, nz, ntimesteps);
  data.speed_of_sound = zeros(nx, ny, nz, ntimesteps);
  data.T = zeros(nx, ny, nz, ntimesteps);
  data.gamma = zeros(nx, ny, nz, ntimesteps);
  data.heat_capacity_at_constant_pressure = zeros(nx, ny, nz, ntimesteps);
  data.species = zeros(nx, ny, nz, ntimesteps, nspecies);

  for k = 1 : ntimesteps
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    timestep = importdata(filename, ' ', 6);
    timepoint_string = timestep.textdata{4};
    time(k) = sscanf(timepoint_string, 't = %f');
    data.rho(:, :, :, k) = reshape(timestep.data(:, 4), nx, ny, nz);
    data.u(:, :, :, k, 1) = reshape(timestep.data(:, 5), nx, ny, nz);
    data.u(:, :, :, k, 2) = reshape(timestep.data(:, 6), nx, ny, nz);
    data.u(:, :, :, k, 3) = reshape(timestep.data(:, 7), nx, ny, nz);
    data.speed_of_sound(:, :, :, k) = reshape(timestep.data(:, 8), nx, ny, nz);
    data.T(:, :, :, k) = reshape(timestep.data(:, 9), nx, ny, nz);
    data.p(:, :, :, k) = reshape(timestep.data(:, 10), nx, ny, nz);
    data.gamma(:, :, :, k) = reshape(timestep.data(:, 11), nx, ny, nz);
    data.heat_capacity_at_constant_pressure(:, :, :, k) = reshape(timestep.data(:, 12), nx, ny, nz);
    data.species(:, :, :, k, :) = reshape(timestep.data(:, 13:end), nx, ny, nz, nspecies);
  end
end
