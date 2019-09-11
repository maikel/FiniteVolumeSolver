function MakeVideo(out, x, time, data)
    v = VideoWriter(out);
    open(v);
    for k = 1:length(time)
       subplot(3, 1, 1);       
       plot(x, data.rho(:, k));
       axis([x(1) x(end) 0 4]);
       ylabel('Density [kg / m^3]');
       title(sprintf('time = %fs', time(k)));
       
       subplot(3, 1, 2);
       plot(x, data.p(:, k)); axis([x(1) x(end) 0 2e5]);
       ylabel('Pressure [Pa]');
       
       subplot(3, 1, 3);
       plot(x, data.T(:, k)); axis([x(1) x(end) 0 3500]);
       ylabel('Temperature [K]');
       xlabel('X-Axis [m]');   
       
       frame = getframe(gcf);
       writeVideo(v, frame);
    end
    close(v);
end

