pathToPlenum = '/group/ag_klima/SFB1029_C01/HLRN/MultiTube/Plenum_x5';
videoPath = '/srv/public/Plenum_x5.avi';

[X, Y, Z, time, data] = ReadPlenumData(pathToPlenum);

%%

f = figure('visible', 'off');
MakeVideo(videoPath, Y, Z, time, data);

function MakeVideo(out, X, Y, time, data)
    v = VideoWriter(out);
    open(v);
    nx = length(X);
    ny = length(Y);
    for k = 1:length(time)
%        subplot(2, 1, 1);       
%        imagesc(X, Y, reshape(data.rho(:, 1, :, k), nx, ny)');
%        axis equal;
%        set(gca, 'Ydir', 'normal');
%        colorbar;
%        caxis([1.0 1.7])
%        ylabel('Density [kg / m^3]');
       fprintf('Write Frame %d.\n', k);
       subplot(2, 1, 1);
       imagesc(X, Y, reshape(data.p(1, :, :, k), nx, ny)');
       axis equal;
       colorbar;
       caxis([8e4 1.3e5])
       set(gca, 'Ydir', 'normal');
       ylabel('Pressure [Pa]');
       title(sprintf('time = %fs', time(k)));
       
       subplot(2, 1, 2);
       imagesc(X, Y, reshape(data.rhou(1, :, :, k, 1) ./ data.rho(1, :, :, k), nx, ny)');
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
