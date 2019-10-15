pathToTube='/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build/IdealGasMix/Matlab';

[X, time, data] = ReadTubeData(pathToTube);

%%
figure(1)
imagesc(X, time * 1000, data.p'); 
set(gca, 'Ydir', 'normal');
xlabel('Position im Rohr [m]');
ylabel('Zeit [ms]');
colorbar;
title('x/t-Diagramm für den Druck [Pa]');

figure(2)
imagesc(X, time, data.rho'); 
set(gca, 'Ydir', 'normal');
xlabel('Position im Rohr [m]');
ylabel('Zeit [ms]');
colorbar;
title('x/t-Diagramm für die Dichte [kg/m3]');

figure(3)
imagesc(X, time, (0.25 * data.species(:, :, 7) ./ (data.species(:, :, 4) / 32))'); 
set(gca, 'Ydir', 'normal');
colorbar;

figure(4)
imagesc(X, time * 1000, (data.u ./ data.speed_of_sound)'); 
set(gca, 'Ydir', 'normal');
xlabel('Position im Rohr [m]');
ylabel('Zeit [ms]');
colorbar;
title('x/t-Diagramm für die Machzahl [-]');

figure(5)
imagesc(X, time, data.T'); 
set(gca, 'Ydir', 'normal');
colorbar;
% title('x/t-Diagramm für Equivalenzverhältnis [-]');
% xaxis('Equivalenzverhältnis H2/O2 [-]');