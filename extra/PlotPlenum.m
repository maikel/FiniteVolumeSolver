pathToPlenum = '/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build_3d/MultiTube/Matlab/Plenum_y0';
videoPath = '/Volumes/Maikel_Intenso/Plenum_y0.avi';

[X, Y, Z, time, data] = ReadPlenumData(pathToPlenum);

%%

f = figure('visible', 'off');
MakeVideo(videoPath, X, Z, time, data);

function MakeVideo(out, X, Y, time, data)
    v = VideoWriter(out);
    open(v);
    nx = length(X);
    ny = length(Y);
    for k = 340:length(time)
      fprintf("k = %d\n", k);
%        subplot(2, 1, 1);       
%        imagesc(X, Y, reshape(data.rho(:, 1, :, k), nx, ny)');
%        axis equal;
%        set(gca, 'Ydir', 'normal');
%        colorbar;
%        caxis([1.0 1.7])
%        ylabel('Density [kg / m^3]');
       
       subplot(2, 1, 1);
       imagesc(X, Y, reshape(data.p(:, 1, :, k), nx, ny)');
       axis equal;
       colorbar;
       caxis([8e4 2e5])
       set(gca, 'Ydir', 'normal');
       ylabel('Pressure [Pa]');
       title(sprintf('time = %fs', time(k)));
       
       subplot(2, 1, 2);
       imagesc(X, Y, reshape(data.rhou(:, 1, :, k, 1) ./ data.rho(:, 1, :, k), nx, ny)');
       axis equal;
       colorbar;
       caxis([0 80])
       set(gca, 'Ydir', 'normal');
       ylabel('Velocity [m/s]');
       xlabel('X-Axis [m]');   

       frame = getframe(gcf);
       writeVideo(v, frame);
    end
    close(v);
end
