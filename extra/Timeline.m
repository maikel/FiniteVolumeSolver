function Timeline(valve, valve_labels)
  figure;
  for j = 1:size(valve, 1)
    for i = 1:size(valve, 2)-1
      patch([valve(j, i), valve(j, i + 1), valve(j, i + 1), valve(j, i)], [j - 1, j - 1, j, j], valve_labels(i));
    end
  end
  axis([0 0.136452 0 5]);
end

