function [X, Y, Z, cycle, time, data] = ReadTimestepData(path)
  fileId = fopen(path);
  size = sscanf(fgetl(fileId), 'size = (%d, %d, %d, %d)');
  dx = sscanf(fgetl(fileId), 'dx = (%f, %f, %f)');
  x0 = sscanf(fgetl(fileId), 'x0 = (%f, %f, %f)');
  time = sscanf(fgetl(fileId), 't = %f');
  cycle = sscanf(fgetl(fileId), 'cycle = %d');
  bin_file = sprintf('%s.bin', path);
  fclose(fileId);
  X = linspace(x0(1), x0(1) + (size(1) - 1)*dx(1), size(1));
  Y = linspace(x0(2), x0(2) + (size(2) - 1)*dx(2), size(2));
  Z = linspace(x0(3), x0(3) + (size(3) - 1)*dx(3), size(3));
  
%   [dir,~,~] = fileparts(path);
  fileId = fopen(bin_file);
  data = reshape(fread(fileId, 'double'), size');
  fclose(fileId);
end

