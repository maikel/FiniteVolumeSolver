pathToPlenums = ["/scratch/guttula/Compressor/Selection/Plenum_x3/", "/scratch/guttula/Compressor/553215/MultiTube_Compressor/Matlab/Plenum_x3/", "/scratch/guttula/Compressor/553402/MultiTube_Compressor/Matlab/Plenum_x3/"];
pattern = '/*.dat';
videoPath = '/scratch/guttula/Compressor/553215/Plenum_x3.avi';

% v = VideoWriter(videoPath);
% open(v);

t_ = [];
maxp_ = [];
meanp_ = [];

for pathToPlenum = pathToPlenums
  fprintf('Process directory %s\n', pathToPlenum);

  files_structure = dir(sprintf('%s%s', pathToPlenum, pattern));
  nFiles = length(files_structure);
  chunkSize = 100;

  t = zeros(nFiles, 1);
  maxp = zeros(nFiles, 1);
  meanp = zeros(nFiles, 1);

  f = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0, 0, 24, 24]);
  for firstFile = 1:chunkSize:nFiles
    [~, Y, Z, time, data] = ReadPlenumDataN(pathToPlenum, firstFile, chunkSize, pattern);
    % MakeVideo(v, Y, Z, time, data, firstFile, nFiles);
    t(firstFile:firstFile+length(time) - 1) = time;
    for k = 1:length(time) 
      maxp(k + firstFile - 1) = max(max(data.p(1, :, :, k)));
      meanp(k + firstFile - 1) = sum(sum(data.p(1, :, :, k) / sum(sum(data.p(1, :, :, k) > 0.0))));
    end
  end

  t_ = [t_; t];
  maxp_ = [maxp_; maxp];
  meanp_ = [meanp_; meanp];
end
% close(v);

f2 = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0, 0, 12, 24]);

subplot(2, 1, 1);
plot(t_, maxp_);

subplot(2, 1, 2);
plot(t_, meanp_);

saveas(f2, '/scratch/guttula/Compressor/553215/Pressure.eps', 'epsc'); 
save('pressure_plots')

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
