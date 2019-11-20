function MakeVideo(out, data)
    v = VideoWriter(out);
    open(v);
    for k = 1:length(data.t)    
       subplot(3, 1, 1);
       imagesc(data.x, data.y, data.rhou(:,:,1,k) ./ data.rho(:,:,k));
       axis equal;
       colorbar;
%        caxis([0 16e5]);
       title(sprintf('Velocity at time = %fs', data.t(k)));
       subplot(3, 1, 2);
       imagesc(data.x, data.y, data.rho(:,:,k));
       axis equal;
       colorbar;
       caxis([0 5]);
       title(sprintf('Density at time = %fs', data.t(k)));
       subplot(3, 1, 3);
       imagesc(data.x, data.y, data.rhoY(:,:,9,k));
       axis equal;
       colorbar;
       caxis([0 5]);
       title(sprintf('Helium at time = %fs', data.t(k)));
       frame = getframe(gcf);
       writeVideo(v, frame);
    end
    close(v);
end

