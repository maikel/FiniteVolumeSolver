function data = ReadHdf5Data_PerfectGas(path)
  raw_data = h5read(path, '/data');  
  sizes = size(raw_data);
  dx = h5readatt(path, '/data', 'cell_size');
  x0 = h5readatt(path, '/data', 'xlower');
  data.x = linspace(x0(1), x0(1) + (sizes(2) - 1)*dx(1), sizes(2));
  data.y = linspace(x0(2), x0(2) + (sizes(1) - 1)*dx(2), sizes(1));
  data.t = h5read(path, '/times');
  nx = length(data.x);
  ny = length(data.y);
  nt = sizes(4);
  data.rho = permute(reshape(raw_data(:,:,1,:), nx, ny, nt), [2 1 3]);
  data.u = permute(reshape(raw_data(:,:,2,:) ./ raw_data(:,:,1,:), nx, ny, nt), [2 1 3]);
  data.v = permute(reshape(raw_data(:,:,3,:) ./ raw_data(:,:,1,:), nx, ny, nt), [2 1 3]);
  data.rhoE = permute(reshape(raw_data(:,:,4,:), nx, ny, nt), [2 1 3]);
  data.p = permute(reshape(raw_data(:,:,5,:), nx, ny, nt), [2 1 3]);
  data.a = permute(reshape(raw_data(:,:,6,:), nx, ny, nt), [2 1 3]);
end
