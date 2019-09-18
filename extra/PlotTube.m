pathToTube='/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build_3d/MultiTube/Matlab/Tube_0';

[X, time, data] = ReadTubeData(pathToTube);

%%
figure(1)
imagesc(X, time, data.p'); 
set(gca, 'Ydir', 'normal');
colorbar;
title('x/t-Diagramm f�r den Druck [Pa]');

figure(2)
imagesc(X, time, data.rho'); 
set(gca, 'Ydir', 'normal');
colorbar;
title('x/t-Diagramm f�r die Dichte [kg/m3]');

figure(3)
imagesc(X, time, (0.25 * data.species(:, :, 7) ./ (data.species(:, :, 4) / 32))'); 
set(gca, 'Ydir', 'normal');
colorbar;
% title('x/t-Diagramm f�r Equivalenzverh�ltnis [-]');
% xaxis('Equivalenzverh�ltnis H2/O2 [-]');