function test_sb_det_simg

clear;
close all;

cur_dir = pwd;

[filename, pathname] = uigetfile({'*.*','All Files (*.*)'},'Import first image');
filename=strcat(char(pathname),char(filename));
img = imread(filename);
cd(cur_dir);

tic;
resolution = 500;
[row,col,dummy] = size(img);
if row > resolution | col > resolution
    if row >= col
        m = resolution/row;
    else
        m = resolution/col;      
    end;
    img = imresize(img,m);    
end
figure,subplot(1,2,1),imshow(img);

[row,col,dummy] = size(img);
img = im2double(img);
map_img = zeros(row,col);
img2 = rgb2ycbcr(img);

% Cb : [.625 .675]
% Cr : [.35 .65]
load sb_color_boundary;

theta = -pi/9;
RotCoeff = [cos(theta) -sin(theta); sin(theta) cos(theta)];

temp1 = img2(:,:,2);
temp2 = img2(:,:,3);
Tdata = [temp1(:) temp2(:)]';
Tdata = RotCoeff*Tdata;

k = 1;
for j=1:col  
    for i=1:row
        if Tdata(1,k)>Cb_low & Tdata(1,k)<Cb_upp & ...
           Tdata(2,k)>Cr_low & Tdata(2,k)<Cr_upp
            map_img(i,j) = 1;
        end    
        k = k+1;
    end
end

subplot(1,2,2),imshow(map_img);
total_time = toc;  
disp(['total time : ',num2str(total_time),'s']);
