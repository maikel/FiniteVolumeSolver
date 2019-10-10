pathToDivider = '/Volumes/Maikel_Intenso/Design5';
videoPath = '/Volumes/Maikel_Intenso/Design5.avi';

[X, Y, time, data] = ReadDividerData(pathToDivider);

%%

f = figure('visible', 'off');
MakeVideo(videoPath, X, Y, time, data);

function MakeVideo(out, X, Y, time, data)
    v = VideoWriter(out);
    open(v);
    for k = 1:length(time)
      fprintf("k = %d\n", k);
     
       imagesc(X, Y, data.p(:, :, k)');
       axis equal;
       colorbar;
       caxis([9e4 6e5])
       set(gca, 'Ydir', 'normal');
       ylabel('Pressure [Pa]');
       title(sprintf('time = %fs', time(k)));
       
       frame = getframe(gcf);
       writeVideo(v, frame);
    end
    close(v);
end
