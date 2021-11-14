function select_sb_color_auto

clear;
close all;

cur_dir = pwd;

cd(['sb color\']);
file_img = dir('*.jpg');
N = length(file_img);
disp(['total images : ',num2str(N)]);

sb_row = 25;
sb_col = 25;

R = []; G = []; B = [];
temp_img1 = [];
temp_img2 = [];
temp_img3 = [];

for i=1:N
    img = imread(file_img(i).name);        
    imshow(img);pause(.5);
    
    img = imresize(img,[sb_row sb_col]);
    img = im2double(img);
    temp1 = img(:,:,1);
    temp2 = img(:,:,2);
    temp3 = img(:,:,3);
    R = [R; temp1(:)];
    G = [G; temp2(:)];
    B = [B; temp3(:)];
    temp_img1 = [temp_img1 temp1];
    temp_img2 = [temp_img2 temp2];
    temp_img3 = [temp_img3 temp3];                
    
    disp(['image: ',num2str(i)]);
end

close all;
sb_dat = [R G B];
sb_dat_img(:,:,1) = temp_img1;
sb_dat_img(:,:,2) = temp_img2;
sb_dat_img(:,:,3) = temp_img3;

cd(cur_dir);
save sb_color_data sb_dat sb_dat_img sb_row sb_col;
disp(['completed..']);
