clear all;
close all;

global params

params.N_Wall = 1000;

%% Geometrie Generation

L_total = 0.15;

p_1 = [0, 0]; 
p_2 = [0, 0.0085]; 
p_3 = [0.03929, 0.00425];
p_4 = [0.03929, 0.0085]; 
% p_5 = [L_total-0.24713, 0.00747]; 
p_6 = [0.08398, 0.01953];
p_7 = [0.08242, 0.02324]; 
p_8 = [0.09602, 0.01953]; 
p_9 = [0.09758, 0.02324]; 
% p_10 = [L_total-0.19987, 0.02430]; 
p_11 = [0.14071, 0.00425];
p_12 = [0.14071, 0.0085];
p_13 = [0.15+0.022, 0];
p_14 = [0.15+0.022, 0.0085];

upper_center_r_up = [0.09, 0.00106]; upper_R_u = 0.02344;
upper_center_r_low = [0.09, 0.00108]; upper_R_l = 0.0194;

upper_div_up_start = p_4(1); upper_div_up_end = p_12(1);
upper_div_low_start = p_3(1); upper_div_low_end = p_11(1);

idx_1_uu = round(params.N_Wall*(p_7(1)-p_4(1))/(upper_div_up_end-upper_div_up_start));
idx_2_uu = round(params.N_Wall*(p_9(1)-p_7(1))/(upper_div_up_end-upper_div_up_start));
idx_3_uu = params.N_Wall-idx_1_uu-idx_2_uu;

idx_1_ul = round(params.N_Wall*(p_6(1)-p_3(1))/(upper_div_low_end-upper_div_low_start));
idx_2_ul = round(params.N_Wall*(p_8(1)-p_6(1))/(upper_div_low_end-upper_div_low_start));
idx_3_ul = params.N_Wall-idx_1_ul-idx_2_ul;

params.upper_div_up(:,1) = [linspace(p_4(1),p_7(1),idx_1_uu) ...
    linspace(p_7(1)+(p_9(1)-p_7(1))/(idx_2_uu+1),p_9(1)-(p_9(1)-p_7(1))/(idx_2_uu+1),idx_2_uu) ...
    linspace(p_9(1),p_12(1),idx_3_uu)];
params.upper_div_up(:,2) = [linspace(p_4(2),p_7(2),idx_1_uu) ...
    linspace(p_7(2)+(p_9(2)-p_7(2))/(idx_2_uu+1),p_9(2)-(p_9(2)-p_7(2))/(idx_2_uu+1),idx_2_uu) ...
    linspace(p_9(2),p_12(2),idx_3_uu)];

idx_start = idx_1_uu+1; idx_end = idx_1_uu+idx_2_uu-1;
params.upper_div_up(idx_start:idx_end,2) = upper_center_r_up(2)+sqrt(upper_R_u^2 -(params.upper_div_up(idx_start:idx_end,1)-upper_center_r_up(1)).^2);

params.upper_div_low(:,1) = [linspace(p_3(1),p_6(1),idx_1_ul) ...
    linspace(p_6(1)+(p_8(1)-p_6(1))/(idx_2_ul+1),p_8(1)-(p_8(1)-p_6(1))/(idx_2_ul+1),idx_2_ul) ...
    linspace(p_8(1),p_11(1),idx_3_ul)];
params.upper_div_low(:,2) = [linspace(p_3(2),p_6(2),idx_1_ul) ...
    linspace(p_6(2)+(p_8(2)-p_6(2))/(idx_2_ul+1),p_8(2)-(p_8(2)-p_6(2))/(idx_2_ul+1),idx_2_ul) ...
    linspace(p_8(2),p_11(2),idx_3_ul)];

idx_start = idx_1_ul+1; idx_end = idx_1_ul+idx_2_ul-1;
params.upper_div_low(idx_start:idx_end,2) = upper_center_r_low(2)+sqrt(upper_R_l^2 -(params.upper_div_low(idx_start:idx_end,1)-upper_center_r_low(1)).^2);


Wall_1(:,1) = linspace(p_1(1),p_13(1),params.N_Wall);
Wall_1(:,2) = linspace(p_1(2),p_13(2),params.N_Wall);

Wall_2(:,1) = linspace(p_3(1),p_11(1),params.N_Wall);
Wall_2(:,2) = linspace(p_3(2),p_11(2),params.N_Wall);

Wall_3(:,1) = params.upper_div_low(:,1); Wall_3(:,2) = params.upper_div_low(:,2);

Wall_4(:,1) = [linspace(p_2(1),p_4(1),params.N_Wall) params.upper_div_up(:,1)' linspace(p_12(1),p_14(1),params.N_Wall)];
Wall_4(:,2) = [linspace(p_2(2),p_4(2),params.N_Wall) params.upper_div_up(:,2)' linspace(p_12(2),p_14(2),params.N_Wall)];



figure(1)

plot(Wall_1(:,1),Wall_1(:,2),'k','Linewidth',2); hold on;
plot(Wall_2(:,1),Wall_2(:,2),'k','Linewidth',2); hold on;
plot(Wall_3(:,1),Wall_3(:,2),'k','Linewidth',2); hold on;
plot(Wall_4(:,1),Wall_4(:,2),'k','Linewidth',2); hold on;

fileID = fopen('wall_1.txt','w'); fprintf(fileID,'%12.8f %12.8f\r\n',Wall_1');
fileID = fopen('wall_2.txt','w'); fprintf(fileID,'%12.8f %12.8f\r\n',Wall_2');
fileID = fopen('wall_3.txt','w'); fprintf(fileID,'%12.8f %12.8f\r\n',Wall_3');
fileID = fopen('wall_4.txt','w'); fprintf(fileID,'%12.8f %12.8f\r\n',Wall_4');
