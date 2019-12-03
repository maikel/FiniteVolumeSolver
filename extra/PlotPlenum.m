pathToPlenums = ["/group/ag_klima/SFB1029_C01/Compressor/Plenum_x1.h5"];
% pattern = '/*.dat';
% videoPath = '/group/ag_klima/SFB1029_C01/Compressor/Plenum_x0.avi';

% v = VideoWriter(videoPath);
% open(v);

t_ = [];
minp_ = [];
maxp_ = [];
meanp_ = [];
variance_ = [];

for pathToPlenum = pathToPlenums
  fprintf('Process HDF5 database %s\n', pathToPlenum);

  info = h5info(pathToPlenum);
  count = info.Datasets(2).Dataspace.Size;
  nFiles = count(end);
  chunkSize = 100;
  chunkCount = count;
  chunkCount(end) = chunkSize;
  first = ones(size(chunkCount));
  X = 1;

  t = zeros(nFiles, 1);
  minp = zeros(nFiles, 1);
  maxp = zeros(nFiles, 1);
  meanp = zeros(nFiles, 1);
  variance = zeros(nFiles, 1);

  f = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0, 0, 24, 24]);
  totalChunks = ceil(nFiles / chunkSize);
  counter = 1;
  for firstFile = 1:chunkSize:nFiles
    fprintf('Process Chunk %d of %d\n', counter, totalChunks);
    counter = counter + 1;
    first(end) = firstFile;
    chunkCount(end) = min(chunkCount(end), nFiles - firstFile + 1);
    [X, ~, ~, time, data] = ReadHdf5DataN(pathToPlenum, first, chunkCount);
    % MakeVideo(v, Y, Z, time, data, firstFile, nFiles);
    t(firstFile:firstFile+length(time) - 1) = time;
    for k = 1:length(time) 
      p = reshape(data.p(:, :, :, 1, k), count(1), count(2));
      mask = p > 0.0;
      minp(k + firstFile - 1) = min(p(mask));
      maxp(k + firstFile - 1) = max(p(mask));
      meanp(k + firstFile - 1) = mean(p(mask));
      variance(k + firstFile - 1) = max(sqrt((p(mask) - mean(p(mask))).^2));
    end
    is = firstFile:firstFile+length(time) - 1;
    stddev = max(sqrt((p - mean(meanp(is)))**2))
    fprintf('stddev %f of mean %f \n', max(variance(firstFile:firstFile+length(time) - 1)), mean(meanp(firstFile:firstFile+length(time) - 1)));
  end

  t_ = [t_; t];
  minp_ = [minp_; minp];
  maxp_ = [maxp_; maxp];
  meanp_ = [meanp_; meanp];
  variance_ = [variance_; variance];
end
% close(v);

f2 = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0, 0, 12, 12]);
hold on;
plot(t_, maxp_);
plot(t_, minp_);
plot(t_, meanp_);
title(sprintf('Druck in x_0 = %fm Entfernung', X(1)));
xlabel('Zeit [s]');
ylabel('Druck [Pa]');
axis([0 0.5 1e5 1.5e5]);

saveas(f2, 'Pressure0.fig');
save('pressure_plots0')

function MakeVideo(v, X, Y, time, data, k0, nFiles)
    nx = length(X);
    ny = length(Y);
    nt = length(time);

    p = reshape(data.p(:,:,:,:), nx, ny, nt);
    level = [0.1e5, linspace(1.2e5,1.3e5, 20)];

    for k = 1:nt
       fprintf('[%3d%%] Write Frame %d of %d.\n', floor(100.0 * (k + k0 - 1) / nFiles), k + k0 - 1, nFiles);

       contourf(X, Y, p(:, :, k)', level);
       c = colorbar;
       c.Label.String = 'Pressure [Pa]';
       caxis([1.2e5 1.3e5]);
       set(gca, 'Ydir', 'normal');
       axis equal;
       ylabel('Z-Axis [m]');
       xlabel('X-Axis [m]');
       title(sprintf('time = %fs', time(k)));

       frame = getframe(gcf);
       writeVideo(v, frame);
    end
end
