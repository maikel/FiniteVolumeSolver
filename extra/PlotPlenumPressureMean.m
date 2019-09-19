pathToPlenum = '/Volumes/Maikel_Intenso/FiniteVolumeSolver_Build_3d/MultiTube/Matlab/Plenum_x0';

[X, Y, Z, time, data] = ReadPlenumData(pathToPlenum);

%%

p = zeros(length(time), 1);

for k = 1:length(time)
  p(k) = sum(sum(data.p(1, :, :, k))) ./ sum(sum(data.p(1, :, :, k) > 0.0));
end

%%
figure
plot(time, p);
axis([time(1) time(end) 3e4 3e5]);
ylabel('Pressure [Pa]');
xlabel('Time [s]');