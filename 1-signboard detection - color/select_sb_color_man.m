function select_sb_color_man

clear;
close all;

cur_dir = pwd;

cd ..
cd('test image\ori\');
file_img = dir('*.jpg');
N = length(file_img);
disp(['total images : ',num2str(N)]);

f1 = figure;
f2 = figure;
 
sb_row = 20;
sb_col = 20;

R = []; G = []; B = [];
temp_img1 = [];
temp_img2 = [];
temp_img3 = [];

for i=1:N
    img = imread(file_img(i).name);        
    figure(f1); imshow(img);
    
    while 1
        figure(f1);
        crop_img = imcrop;
        figure(2); imshow(crop_img);
        Ask = questdlg('Select This Region?','Genie Question','Yes','No','Cancel','Yes');
        switch Ask
            case 'Yes'                
                crop_img = imresize(crop_img,[sb_row sb_col]);
                crop_img = im2double(crop_img);
                temp1 = crop_img(:,:,1);
                temp2 = crop_img(:,:,2);
                temp3 = crop_img(:,:,3);
                R = [R; temp1(:)];
                G = [G; temp2(:)];
                B = [B; temp3(:)];
                temp_img1 = [temp_img1 temp1];
                temp_img2 = [temp_img2 temp2];
                temp_img3 = [temp_img3 temp3];                
                figure(2); clf;
            case 'No'
                figure(2); clf;
            otherwise
                figure(2); clf;
                break;
        end
    end
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
