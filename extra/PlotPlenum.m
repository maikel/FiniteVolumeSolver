pathToPlenum = '/Volumes/Maikel_Intenso/Selection/Plenum_x6';
pattern = '/*.dat';
videoPath = '/Users/maikel/Plenum.avi';

files_structure = dir(sprintf('%s%s', pathToPlenum, pattern));
nFiles = 100;
chunkSize = 100;

f = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0, 0, 24, 24]);
v = VideoWriter(videoPath);
open(v);
for firstFile = 1:chunkSize:nFiles
  [~, Y, Z, time, data] = ReadPlenumDataN(pathToPlenum, firstFile, chunkSize, pattern);
  MakeVideo(v, Y, Z, time, data);
end
close(v);

function MakeVideo(v, X, Y, time, data)
    nx = length(X);
    ny = length(Y);
    nt = length(time);

    u = reshape(data.rhou(:,:,:,:,1) ./ data.rho, nx, ny, nt);
    a = reshape(data.speed_of_sound, nx, ny, nt);

    Mach = u ./ a;

    levels = [-1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

    for k = 1:nt
       fprintf('Write Frame %d.\n', k);

       contourf(X, Y, Mach(:, :, k)', levels);
       set(gca, 'Ydir', 'normal');
       axis equal;
       c = colorbar;
       c.Label.String = 'Mach number [-]';
       caxis([0.0 1.0])
       ylabel('Z-Axis [m]');
       xlabel('X-Axis [m]');
       title(sprintf('time = %fs', time(k)));

       frame = getframe(gcf);
       writeVideo(v, frame);
    end
end
