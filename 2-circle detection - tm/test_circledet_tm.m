function test_circledet_tm

clear;
close all;
clc;

[filename, pathname] = uigetfile({'*.*','All Files (*.*)'},'Import first image');
filename=strcat(char(pathname),char(filename));
I=imread(filename);

% I = imread('blood1.tif');
% I = imread('eight.tif');

resolution = 500;
[row_ori,col_ori,dummy] = size(I);
if row_ori > resolution | col_ori > resolution
    if row_ori >= col_ori
        m = resolution/row_ori;
    else
        m = resolution/col_ori;      
    end;
    I = imresize(I,m);    
    [row_ori,col_ori,dummy] = size(I);
end

if isrgb(I)
    I = rgb2gray(I);
end
figure,subplot(1,2,1),imshow(I);
title(['original image']);
drawnow;

% BW = im2bw(I,graythresh(I));
BW = ~im2bw(I,graythresh(I));
subplot(1,2,2),imshow(BW); title(['thresholding']); pause(1);
SE = strel('square',5);
BW = imdilate(BW,SE);
subplot(1,2,2),imshow(BW); title(['dilate']); pause(1);
BW = imfill(BW, 'holes');
subplot(1,2,2),imshow(BW); title(['fill holes']); pause(1);
BW = imerode(BW,SE);
subplot(1,2,2),imshow(BW); title(['eroding']); pause(1);

tic;
%-------------------------------------------%
load circle_filter;
box_color = ['green'];
center_color = ['yellow'];
num_scales = size(fft_cir,1);
fft_res = size(fft_cir{1},1);
cir_threshold = .8;
max_corr = 232.6;  
%-------------------------------------------%

cand_img = BW;
m = 1;
if row_ori > fft_res | col_ori > fft_res
    if row_ori >= col_ori
        m = fft_res/row_ori;
    else
        m = fft_res/col_ori;      
    end;
    cand_img = imresize(cand_img,m,'bicubic');
end

cand_img = im2double(cand_img);
[row,col] = size(cand_img);

srow = row;
scol = col;

fft_cand_img = fft2(cand_img,fft_res,fft_res);
for scale = 1:num_scales
    res = scale_ratio^(scale-1);
    if scale == 1
        fft_cur = fft_cand_img;
    else
        srow = round(row/res);
        scol = round(col/res);
        fft_cur = fnn_downs_fft(fft_cand_img,res);
    end 
    if (srow > cir_height) & (scol > cir_width)
        fft_scale_cir = fft_cir{scale};        
        [output_val,output_pos,output_box] = locate_cir(row_ori,col_ori,m,srow,scol,...
                                                        fft_cur,fft_scale_cir,...
                                                        res,cir_height,cir_width,...
                                                        cir_threshold,max_corr);
        m_val(scale,1) = {output_val};
        m_pos(scale,1) = {output_pos};
        m_box(scale,1) = {output_box};     
    else
        break;
    end
    disp(['searching on scale: ',num2str(scale),' completed']);
end

figure,imshow(I);
[output_val,output_pos,output_box,output_agg] = post_processing(m_val,m_pos,m_box);
draw_box(output_val,output_pos,output_box,output_agg,box_color,center_color);
total_time = toc;  
disp(['total time : ',num2str(total_time),'s']);

%############################################################################################
function [output_val,output_pos,output_box] = locate_cir(row_ori,col_ori,m,trow,tcol,...
                                                         fft_cand,fft_cir,...
                                                         res,cir_height,cir_width,...
                                                         cir_threshold,max_corr)

[fft_res,dummy] = size(fft_cir);
output = real(ifft2(fft_cand.*fft_cir));
output = output(cir_height:trow,cir_width:tcol);
map = -1*ones(trow,tcol);
map(cir_height:trow,cir_width:tcol) = output/max_corr;

[i,j] = find(map>=cir_threshold);
output_val = [];
output_pos = [];
output_box = [];
num_detect = length(i);
thet = cos(pi/4);
for k=1:num_detect
    row = round(i(k)*res/m);
    col = round(j(k)*res/m);
    ColumnStart = max(1,col-(cir_width*res/m));
    RowStart = max(1,row-(cir_height*res/m));
    ColumnEnd = min(col_ori,col);
    RowEnd = min(row_ori,row);    
    center_row = round((RowEnd-RowStart)/2)+RowStart;
    center_col = round((ColumnEnd-ColumnStart)/2)+ColumnStart;

    pyth_y = round((ColumnEnd-ColumnStart)*thet/2);
    RowStart = max([1 (center_row-pyth_y)]);
    RowEnd = min([row_ori (center_row+pyth_y)]);
    ColumnStart = max([1 (center_col-pyth_y)]);
    ColumnEnd = min([col_ori (center_col+pyth_y)]);

    output_val = [output_val; map(i(k),j(k))];
    output_pos = [output_pos; center_row center_col];
    output_box = [output_box; RowStart RowEnd ColumnStart ColumnEnd];
end

%############################################################################################
function draw_box(out_val,out_pos,out_box,out_agg,box_color,center_color)

num_detect = length(out_val);
if num_detect > 0
    for k=1:num_detect
        RowStart = out_box(k,1);
        RowEnd = out_box(k,2);
        ColumnStart = out_box(k,3);
        ColumnEnd = out_box(k,4);
        
        set(line([ColumnStart, ColumnEnd], [RowStart, RowStart]), 'Color', box_color, 'Linewidth', 1);        
        set(line([ColumnStart, ColumnEnd], [RowEnd, RowEnd]), 'Color', box_color, 'Linewidth', 1);        
        set(line([ColumnStart, ColumnStart], [RowStart, RowEnd]), 'Color', box_color, 'Linewidth', 1);        
        set(line([ColumnEnd, ColumnEnd], [RowStart, RowEnd]), 'Color', box_color, 'Linewidth', 1);      
        
        center_row = out_pos(k,1);
        center_col = out_pos(k,2);
        
        set(line([center_col, center_col],[(center_row-2), (center_row+2)]),'color',center_color,'linewidth',1);
        set(line([(center_col-2), (center_col+2)],[center_row, center_row]),'color',center_color,'linewidth',1);
        drawnow;
        disp(['score : ',num2str(out_val(k)),' aggregate : ',num2str(out_agg(k))]);        
    end  
else
    disp(['No Circle Detected!']);
end

%############################################################################################
function [out_val,out_pos,out_box,out_agg] = post_processing(m_val,m_pos,m_box)

%------------------------------------------------------%
m_dis2d = 5;
m_dis3d = 5;
thres_2d = 2;
thres_3d = 2;
%------------------------------------------------------%

c_val = [];
c_pos = [];
c_box = [];
c_agg = [];
out_val = [];
out_pos = [];
out_box = [];
out_agg = [];

num_scale = length(m_val);
for scale=1:num_scale    
    s_val = m_val{scale};
    num_detect = length(s_val);
    
    if num_detect ~= 0;
        s_pos = m_pos{scale};   
        s_box = m_box{scale}; 
        [dummy,ind] = sort(s_val);
        ind = flipud(ind);
        mark = zeros(num_detect,1);
        counter = 0;
        mark_cand = 1;
        
        while counter < num_detect
            cand_highest = 0;
            for k=1:num_detect
                if mark(k) == 0 
                    if cand_highest == 0
                        mark(k) = mark_cand;
                        mark_cand = mark_cand+1;
                        cand_highest = 1;
                        counter = counter+1;
                        pos = s_pos(ind(k),:);
                        m_dis2d = (abs(pos(1)-s_box(ind(k),1))+abs(pos(2)-s_box(ind(k),3)))/2;
                        m_dis2d = 2*m_dis2d/3;
                    else 
                        multi_dist = abs(pos-s_pos(ind(k),:));
                        if (multi_dist(1) <= m_dis2d) & (multi_dist(2) <= m_dis2d)    
                            mark(k) = mark_cand-1;
                            counter = counter+1;
                        end
                    end
                end
            end
        end
        
        for k=1:mark_cand-1
            [val,cand_ind] = find(mark == k);
            aggregate = length(val);
            %if aggregate >= thres_2d
                for i=1:num_detect
                    if mark(i) == k
                        c_val = [c_val; s_val(ind(i))];
                        c_pos = [c_pos; s_pos(ind(i),:)];
                        c_box = [c_box; s_box(ind(i),:)];
                        c_agg = [c_agg; aggregate];
                        break;
                    end
                end
            %end
        end    
    end
end

num_detect = length(c_val);

if num_detect ~= 0;
    [dummy,ind] = sort(c_val);
    ind = flipud(ind);
    mark = zeros(num_detect,1);
    counter = 0;
    mark_cand = 1;
    
    while counter < num_detect
        cand_highest = 0;
        for k=1:num_detect
            if mark(k) == 0 
                if cand_highest == 0
                    mark(k) = mark_cand;
                    mark_cand = mark_cand+1;
                    cand_highest = 1;
                    counter = counter+1;
                    pos = c_pos(ind(k),:);
                    m_dis3d = (abs(pos(1)-c_box(ind(k),1))+abs(pos(2)-c_box(ind(k),3)))/2;
                    m_dis3d = 2*m_dis3d/3;
                else 
                    multi_dist = abs(pos-c_pos(ind(k),:));
                    if (multi_dist(1) <= m_dis3d) & (multi_dist(2) <= m_dis3d)
                        mark(k) = mark_cand-1;
                        counter = counter+1;
                    end
                end
            end
        end
    end
    
    for k=1:mark_cand-1
        cand_highest = 0;
        acc_agg = 0;
        for i=1:num_detect
            if mark(i) == k
                if cand_highest == 0
                    cand_highest = 1;
                    acc_agg = c_agg(ind(i));
                    temp_val = c_val(ind(i));
                    temp_pos = c_pos(ind(i),:);
                    temp_box = c_box(ind(i),:);
                else
                    acc_agg = acc_agg + c_agg(ind(i));
                end       
            end
        end
        
        if (cand_highest ~= 0) & acc_agg >= thres_3d
            out_val = [out_val; temp_val];
            out_pos = [out_pos; temp_pos];
            out_box = [out_box; temp_box];
            out_agg = [out_agg; acc_agg];
        end            
    end    
end

%############################################################################################
function fft_ds_img = fnn_downs_fft(fft_img,res)

[r,c] = size(fft_img);

row_range = round(r/res);
col_range = round(c/res);

row_s = round((r-row_range)/2)+1;
row_e = row_s+row_range-1;
col_s = round((c-col_range)/2)+1;
col_e = col_s+col_range-1;

F = fftshift(fft_img);
fft_ds_img = F(row_s:row_e,col_s:col_e)/(res^2);
fft_ds_img = ifftshift(fft_ds_img);