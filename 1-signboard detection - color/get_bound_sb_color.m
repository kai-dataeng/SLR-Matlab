function get_bound_sb_color

clc;
close all;
clear;
load sb_color_data;

percent_std_x = 80;
percent_std_y = 70;
str = ['b.']; 

temp = rgb2ycbcr(sb_dat_img);
temp1 = temp(:,:,1);
temp2 = temp(:,:,2);
temp3 = temp(:,:,3);
sb_ycbcr = [temp1(:) temp2(:) temp3(:)];

figure,plot(sb_ycbcr(:,2),sb_ycbcr(:,3),'r.');
axis([0 1 0 1]);
title('CbCr space');
xlabel('Cb'); 
ylabel('Cr');
grid on;

Pdata = sb_ycbcr(:,2:3)';
[d,M] = size(Pdata);

theta = -pi/9;
RotCoeff = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Pdata = RotCoeff*Pdata;

avg_x = mean(Pdata(1,:));
std_x = std(Pdata(1,:));
avg_y = mean(Pdata(2,:));
std_y = std(Pdata(2,:));

ct = std_x*6/2;
mn1 = avg_x-ct; mx1 = avg_x+ct;
mm = minmax([mn1 mx1]);
res = (mm(2)-mm(1))/100;
x = [mm(1):res:mm(2)];

ct = std_y*6/2;
mn2 = avg_y-ct; mx2 = avg_y+ct;
mm = minmax([mn2 mx2]);
res = (mm(2)-mm(1))/100;
y = [mm(1):res:mm(2)];

p_x = normpdf(x,avg_x,std_x);
p_y = normpdf(y,avg_y,std_y);

figure,plot(Pdata(1,:),Pdata(2,:),'r.');
axis([0 1 0 1]);
title('CbCr space (AFTER ROTATION)');
xlabel('Cb'); 
ylabel('Cr');
grid on;

hold on;
A_ND = 0.25;
p_x = A_ND*p_x/max(p_x);
area(x,p_x,'FaceColor','yellow','EdgeColor',str(1),'LineWidth',2);

p_y = A_ND*p_y/max(p_y);
area(p_y,y,'FaceColor','yellow','EdgeColor',str(1),'LineWidth',2);

set(line([avg_x avg_x],[0 1]),'Color',[0 0 0]);
set(line([0 1],[avg_y avg_y]),'Color',[0 0 0]);

h_min = -3*percent_std_x*std_x/100+avg_x;
h_max = 3*percent_std_x*std_x/100+avg_x;

set(line([h_min h_min],[0 1]),'Color','green','LineWidth',2);
set(line([h_max h_max],[0 1]),'Color','green','LineWidth',2);

v_min = -3*percent_std_y*std_y/100+avg_y;
v_max = 3*percent_std_y*std_y/100+avg_y;

% set(line([0 1],[v_min v_min]),'Color','green','LineWidth',2);
% set(line([0 1],[v_max v_max]),'Color','green','LineWidth',2);

hold off;
drawnow;

display([' Cb: [',num2str(h_min),' ',num2str(h_max),']']);
display([' Cr: [',num2str(v_min),' ',num2str(v_max),']']);

% Cb_low = h_min;
% Cb_upp = h_max;
% Cr_low = v_min;
% Cr_upp = v_max;

% Cb : [.625 .675]
% Cr : [.35 .65]
Cb_low = .62;
Cb_upp = .69;
Cr_low = .32;
Cr_upp = .65;

save sb_color_boundary Cb_low Cb_upp Cr_low Cr_upp