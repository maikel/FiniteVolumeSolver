function [X, time, data, columns] = ReadTubeData(directory_path)
  files_structure = dir(sprintf('%s%s', directory_path, '/*.dat'));
  if (length(files_structure) == 0) 
     ME = MException('ReadTubeData:noFilesFound', ...
        'Could not find any *.dat files in directory %s.', directory_path);
     throw(ME);
  end
  all_filenames = {files_structure(:).name};
  timestep = importdata(sprintf('%s/%s', directory_path, all_filenames{1}), ' ', 4);
  nx = size(timestep.data, 1);

  columns = timestep.colheaders;
  X = timestep.data(:, 1);
  time = zeros(length(all_filenames), 1);
  data.rho = zeros(nx, length(all_filenames));
  data.u = zeros(nx, length(all_filenames));
  data.p = zeros(nx, length(all_filenames));
  data.T = zeros(nx, length(all_filenames));

  for k = 1 : length(all_filenames)
    filename = sprintf('%s/%s', directory_path, all_filenames{k});
    timestep = importdata(filename, ' ', 4);
    timepoint_string = timestep.textdata{2};
    time(k) = sscanf(timepoint_string, 't = %f');
    data.rho(:, k) = timestep.data(:, 2);
    data.u(:, k) = timestep.data(:, 3);
    data.T(:, k) = timestep.data(:, 4);
    data.p(:, k) = timestep.data(:, 5);
  end
end

