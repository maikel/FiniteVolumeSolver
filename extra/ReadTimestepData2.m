function [X, Y, cycle, time, data] = ReadTimestepData2(path)
  fileId = fopen(path);
  size = sscanf(fgetl(fileId), 'size = (%d, %d, %d)');
  dx = sscanf(fgetl(fileId), 'dx = (%f, %f)');
  x0 = sscanf(fgetl(fileId), 'x0 = (%f, %f)');
  time = sscanf(fgetl(fileId), 't = %f');
  cycle = sscanf(fgetl(fileId), 'cycle = %d');
  bin_file = sscanf(fgetl(fileId), 'data_file = %s');
  fclose(fileId);
  X = linspace(x0(1), x0(1) + (size(1) - 1)*dx(1), size(1));
  Y = linspace(x0(2), x0(2) + (size(2) - 1)*dx(2), size(2));
  
  [dir,~,~] = fileparts(path);
  fileId = fopen(sprintf('%s/%s', dir, bin_file));
  data = reshape(fread(fileId, 'double'), size');
  fclose(fileId);
end

