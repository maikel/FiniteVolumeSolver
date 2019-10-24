pathToPlenum = '/Volumes/Maikel_Intenso/Selection_No/Plenum_x6';

[X, Y, Z, time, data] = ReadPlenumDataN(pathToPlenum, 1000, 2000);

%%

p = zeros(length(time), 1);

for k = 1:length(time)
  p(k) = max(max(data.p(1, :, :, k))) / 1.0e5;
end

%%
% figure('visible', 'off');
t_ms = time * 1000;
f = figure('Units', 'centimeters', 'Position', [0, 0, 11.5, 6]);
plot(t_ms, p);
axis([10.0 30.0 0.9 4.0]);
ylabel('Druck [Bar]', 'FontSize', 12);
xlabel('Zeit [ms]', 'FontSize', 12);
% saveas(gcf,'PressureOverTime.png')

fprintf('Max Pressure: %f Bar\n', max(p));
