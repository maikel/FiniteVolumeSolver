function data = ReadHdf5Data_IdealGas(path)
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
  nvars = sizes(3);
  ndim = 3 - (21 - nvars); 
  irhoE = 2+ndim;
  iP = irhoE + 1 + 11;
  data.rho = permute(reshape(raw_data(:,:,1,:), nx, ny, nt), [2 1 3]);
  data.rhou = permute(reshape(raw_data(:,:,2:irhoE-1,:), nx, ny, ndim, nt), [2 1 3 4]);
  data.rhoE = permute(reshape(raw_data(:,:,irhoE,:), nx, ny, nt), [2 1 3]);
  data.rhoY = permute(reshape(raw_data(:,:,irhoE+1:iP-1,:), nx, ny, 11, nt), [2 1 3 4]);
  data.p = permute(reshape(raw_data(:,:,iP,:), nx, ny, nt), [2 1 3]);
  data.T = permute(reshape(raw_data(:,:,iP+1,:), nx, ny, nt), [2 1 3]);
  data.a = permute(reshape(raw_data(:,:,iP+2,:), nx, ny, nt), [2 1 3]);
  data.c_p = permute(reshape(raw_data(:,:,iP+3,:), nx, ny, nt), [2 1 3]);
  data.gamma = permute(reshape(raw_data(:,:,iP+4,:), nx, ny, nt), [2 1 3]);

end
