pathToPlenum = '/group/ag_klima/SFB1029_C01/HLRN/MultiTube/Plenum_x5/';

[X, Y, Z, time, data] = ReadPlenumData(pathToPlenum);

%%

p = zeros(length(time), 1);

for k = 1:length(time)
  p(k) = sum(sum(data.p(1, :, :, k))) ./ sum(sum(data.p(1, :, :, k) > 0.0));
end

%%
figure('visible', 'off');
plot(time, p);
axis([time(1) time(end) 3e4 3e5]);
ylabel('Pressure [Pa]');
xlabel('Time [s]');
saveas(gcf,'PressureOverTime.png')

fprintf('Mean Pressure: %f Pa\n', mean(p));
