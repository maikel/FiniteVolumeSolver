function [X, Y, Z, time, data] = ReadPlenumData(directory_path, pattern)
  if nargin == 1
    pattern = '/*.dat';
  end
  files_structure = dir(sprintf('%s%s', directory_path, pattern));
  ntimesteps = length(files_structure);
  if (ntimesteps == 0) 
     ME = MException('ReadPlenumData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  filename = sprintf('%s/%s', directory_path, all_filenames{1});
  [X, Y, Z, ~, ~, timestep] = ReadTimestepData(filename);
  nspecies = size(timestep, 4) - 10;
    
  nx = length(X);
  ny = length(Y);
  nz = length(Z);
  
  time = zeros(ntimesteps, 1);
  data.rho = zeros(nx, ny, nz, ntimesteps);
  data.rhou = zeros(nx, ny, nz, ntimesteps, 3);
  data.rhoE = zeros(nx, ny, nz, ntimesteps);
  data.speed_of_sound = zeros(nx, ny, nz, ntimesteps);
  data.T = zeros(nx, ny, nz, ntimesteps);
  data.p = zeros(nx, ny, nz, ntimesteps);
  data.gamma = zeros(nx, ny, nz, ntimesteps);
  data.heat_capacity_at_constant_pressure = zeros(nx, ny, nz, ntimesteps);
  data.species = zeros(nx, ny, nz, ntimesteps, nspecies);

  for k = 1 : ntimesteps
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    [~, ~, ~, ~, t, timestep] = ReadTimestepData(filename);
    time(k) = t;
    
    data.rho(:, :, :, k) = timestep(:, : , :, 1);
    data.rhou(:, :, :, k, :) = timestep(:, : , :, 2:4);
    data.rhoE(:, :, :, k) = timestep(:, : , :, 5);
   
    sEnd = 5 + nspecies;
    data.species(:, :, :, k, :) = timestep(:, : , :, 6:sEnd);
    
    data.p(:, :, :, k) = timestep(:, : , :, sEnd+1);
    data.T(:, :, :, k) = timestep(:, : , :, sEnd+2);
    data.speed_of_sound(:, :, :, k) = timestep(:, : , :, sEnd+3);
    data.heat_capacity_at_constant_pressure(:, :, :, k) = timestep(:, : , :, sEnd+4);
    data.gamma(:, :, :, k) = timestep(:, : , :, sEnd+5);
  end
end
