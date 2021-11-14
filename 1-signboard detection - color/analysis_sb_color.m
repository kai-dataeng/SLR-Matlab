function analysis_sb_color
close all;
clear;
load sb_color_data;

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

figure,plot(sb_ycbcr(:,1),sb_ycbcr(:,2),'r.');
axis([0 1 0 1]);
title('YCb space');
xlabel('Y'); 
ylabel('Cb');
grid on;

figure,plot(sb_ycbcr(:,1),sb_ycbcr(:,3),'r.');
axis([0 1 0 1]);
title('YCr space');
xlabel('Y'); 
ylabel('Cr');
grid on;

figure,plot3(sb_dat(:,1),sb_dat(:,2),sb_dat(:,3),'r.');
axis([0 1 0 1 0 1]);
title('RGB space');
xlabel('R'); 
ylabel('G');
zlabel('B');
grid on;