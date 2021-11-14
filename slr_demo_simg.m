function slr_demo_simg

warning off;
clc;
clear;
close all;

cur_dir = pwd;
load('D:\WORK\w-MATLAB\speed limit recognition\1-signboard detection - color\sb_color_boundary');
load('D:\WORK\w-MATLAB\speed limit recognition\2-circle detection - tm\circle_filter');
load('D:\WORK\w-MATLAB\speed limit recognition\3-digit recognition - nn\slr_nn');
load slr_init_data;

% cd(['test image\']);
cd('speed-limit - dataset\image sequence\');
[filename, pathname] = uigetfile({'*.*','All Files (*.*)'},'Import first image');
file_img = strcat(char(pathname),char(filename));

%++++++++++++++++++++++++++++++++++++%
disp_time = 1;
%++++++++++++++++++++++++++++++++++++%

tic;
result = 0;
img = imread(file_img);

[row_ori,col_ori,dummy] = size(img);

if(row_ori > min_res && col_ori > min_res)
    h_upp = round(h_end*col_ori/h_bin);
    h_low = round(h_sta*col_ori/h_bin);
    img(:,h_low:h_upp,:) = [];
    v_upp = round(v_end*row_ori/v_bin);
    v_low = round(v_sta*row_ori/v_bin);
    img(1:v_low,:,:) = [];
    img(v_upp:end,:,:) = [];
end

if row_ori > resolution || col_ori > resolution
    if row_ori >= col_ori
        m = resolution/row_ori;
    else
        m = resolution/col_ori;
    end
    img = imresize(img,m);
    [row_ori,col_ori,dummy] = size(img);
end

figure,imshow(img);
title('original image');
drawnow;
pause(disp_time);

%------------- Sign board border detection -------------%
map_img = zeros(row_ori,col_ori);
img = im2double(img);
img2 = rgb2ycbcr(img);
RotCoeff = [cos(theta) -sin(theta); sin(theta) cos(theta)];

temp1 = img2(:,:,2);
temp2 = img2(:,:,3);
Tdata = [temp1(:) temp2(:)]';
Tdata = RotCoeff*Tdata;

k = 1;
for j=1:col_ori
    for i=1:row_ori
        if Tdata(1,k)>Cb_low && Tdata(1,k)<Cb_upp && ...
                Tdata(2,k)>Cr_low && Tdata(2,k)<Cr_upp
            map_img(i,j) = 1;
        end
        k = k+1;
    end
end

figure,imshow(map_img);
title('sign board border detection');
drawnow;
pause(disp_time);

%------------- Circle detection -------------%
BW = logical(map_img);

figure,imshow(BW);
title('circle detection');
drawnow;
pause(disp_time);

cand_img = BW;
m = 1;
if row_ori > fft_res || col_ori > fft_res
    if row_ori >= col_ori
        m = fft_res/row_ori;
    else
        m = fft_res/col_ori;
    end
    cand_img = imresize(cand_img,m,'bicubic');
end

cand_img = im2double(cand_img);
[row,col] = size(cand_img);

srow = row;
scol = col;

m_val = {}; m_pos = {}; m_box = {};
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
    if (srow > cir_height) && (scol > cir_width)
        fft_scale_cir = fft_cir{scale};
        [output_val,output_pos,output_box] = locate_cir(row_ori,col_ori,m,srow,scol,...
            fft_cur,fft_scale_cir,...
            res,cir_height,cir_width,...
            cir_threshold,max_corr);
        m_val(scale,1) = {output_val};
        m_pos(scale,1) = {output_pos};
        m_box(scale,1) = {output_box};
    else
        break
    end
    %             disp(['searching on scale: ',num2str(scale),' completed']);
end

SE = strel('square',5);
BW = imdilate(BW,SE);
BW = imfill(BW, 'holes');
BW = imerode(BW,SE);

[output_val,output_pos,output_box,output_agg] = post_processing(m_val,m_pos,m_box);
draw_box(output_val,output_pos,output_box,output_agg,box_color,center_color);

%------------- Estimate and crop digits area -------------%
num_detect = length(output_val);
if num_detect > 0
    pause(disp_time);
    RowStart = output_box(1,1);
    RowHeight = output_box(1,2)-RowStart;
    ColumnStart = output_box(1,3);
    ColumnWidth = output_box(1,4)-ColumnStart;
    
    if ndims(img) == 3
        img = rgb2gray(img);
    end
    
    img = imcrop(img,[ColumnStart RowStart ColumnWidth RowHeight]);
    da_img = imresize(img,[da_size da_size],'bicubic');
    level = graythresh(da_img);
    BW = ~im2bw(da_img,level);
    figure,imshow(BW);
    title('Estimate digits area');
    drawnow;
    pause(disp_time);
    
    BW = imclearborder(BW);
    imshow(BW);
    title('remove unwanted border');
    drawnow;
    pause(disp_time);
    
    BW = bwareaopen(BW,uw_area);
    imshow(BW);
    title('remove unwanted objects');
    drawnow;
    pause(disp_time);
    
    BW_label = bwlabel(BW);
    num_obj = max(BW_label(:));
    
    temp_map = temp_map*0;
    obj_ind = []; obj_area = []; obj_ord = []; digit_data = [];
    
    if num_obj<2
        imshow([]);
        title('object less than 2 objects');
        drawnow;
    else
        for k=1:num_obj
            [r,c] = find(BW_label == k);
            obj_area = [obj_area; length(r)];
            obj_ind = [obj_ind; struct('row_ind',r,'col_ind',c)];
        end
        
        min_area = .3*max(obj_area);
        obj_count = 0;
        for k=1:num_obj
            if (obj_area(k)>min_area)
                for i=1:obj_area(k)
                    temp_map(obj_ind(k).row_ind(i),obj_ind(k).col_ind(i)) = 1;
                end
                
                columnWeight = col_mat.*temp_map;
                cMean = round(sum(columnWeight(:))/obj_area(k));
                a1 = sum(columnWeight(:).^2);
                a2 = (sum(columnWeight(:))^2)/obj_area(k);
                cStd = round(sqrt((a1-a2)/(obj_area(k)-1)));
                rowWeight = row_mat.*temp_map;
                rMean = round(sum(rowWeight(:))/obj_area(k));
                a1 = sum(rowWeight(:).^2);
                a2 = (sum(rowWeight(:))^2)/obj_area(k);
                rStd = round(sqrt((a1-a2)/(obj_area(k)-1)));
                rowStd = rStd*std_weight;
                colStd = cStd*std_weight;
                
                ColumnStart = max(1, cMean - min(colStd, ceil(da_size/2)));
                ColumnEnd = min(da_size, cMean + min(colStd, ceil(da_size/2)));
                RowStart = max(1, rMean - min(rowStd, ceil(da_size/2)));
                RowEnd = min(da_size, rMean + min(rowStd, ceil(da_size/2)));
                
                figure,imshow(temp_map); title('detected digit');
                set(line([ColumnStart, ColumnEnd], [RowStart, RowStart]), 'Color', 'red', 'Linewidth', 2);
                set(line([ColumnStart, ColumnEnd], [RowEnd, RowEnd]), 'Color', 'red', 'Linewidth', 2);
                set(line([ColumnStart, ColumnStart], [RowStart, RowEnd]), 'Color', 'red', 'Linewidth', 2);
                set(line([ColumnEnd, ColumnEnd], [RowStart, RowEnd]), 'Color', 'red', 'Linewidth', 2);
                drawnow;
                pause(disp_time);
                
                ColumnWidth = max(1,(ColumnEnd-ColumnStart));
                RowHeight = max(1,(RowEnd-RowStart));
                obj_hscale = RowHeight/ColumnWidth;
                
                if obj_hscale>min_obj_hscale & obj_count<4
                    
                    obj_count = obj_count+1;
                    obj_ord = [obj_ord; ColumnStart];
                    
                    if obj_hscale>=3
                        rStd = rStd*(std_weight+.1);
                        cStd = cStd*(std_weight+10);
                        ColumnStart = max(1, cMean - min(cStd, ceil(da_size/2)));
                        ColumnEnd = min(da_size, cMean + min(cStd, ceil(da_size/2)));
                        RowStart = max(1, rMean - min(rStd, ceil(da_size/2)));
                        RowEnd = min(da_size, rMean + min(rStd, ceil(da_size/2)));
                        ColumnWidth = max(1,(ColumnEnd-ColumnStart));
                        RowHeight = max(1,(RowEnd-RowStart));
                        digit = imcrop(temp_map,[ColumnStart RowStart ColumnWidth RowHeight]);
                    else
                        digit = imcrop(temp_map,[ColumnStart RowStart ColumnWidth RowHeight]);
                        digit = remove_empty_background(digit);
                    end
                    digit = imresize(digit,[slr_height slr_width]);
                    digit_data = [digit_data; struct('data',digit)];
                end
                temp_map = temp_map*0;
            end
        end
        
        if obj_count>1
            [dummy,ind] = sort(obj_ord);
            temp = []; P = [];
            for i=1:length(ind)
                xx =  digit_data(ind(i)).data;
                temp = [temp xx];
                P = [P xx(:)];
            end
            digit_img = .5*ones(slr_height+20,(obj_count*slr_width+20));
            digit_img(11:(10+slr_height),11:(obj_count*slr_width+10)) = temp;
            figure,imshow(digit_img);
            title('digit arrangement');
            drawnow;
            
            [result,str] = digit_recog(net,P,nn_threshold);
        end
    end
else
    imshow([]);
    title('No Sign Board Detected');
    drawnow;
end

total_time = toc;
disp(['total time : ',num2str(total_time),'s']);

if result
    disp(['RESULT for image ',file_img,' : ',str,' km/h']);
else
    disp(['RESULT for image ',file_img,' : ------------']);
end

disp(['completed...']);
cd(cur_dir);

%############################################################################################
function [output,str] = digit_recog(net,data,threshold);

N = size(data,2);
Y = sim(net,data);

if size(Y,1)>10
    Y = Y(1:(end-1),:);
end

str = ['']; 
output = 1;
for i=1:N    
    [max_score,Pind] = max(Y(:,i));
    if max_score >= threshold
        str = strcat(str,[num2str(Pind-1)]);
    else
        str = ['unrecognized digit!!'];
        output = 0;
        break;
    end        
end

%############################################################################################
function digit = remove_empty_background(digit)

dist_empty = 2;

cc = sum(digit,1);
rr = sum(digit,2);

ColumnStart = 1;  ColumnEnd = size(digit,2);
RowStart = 1;  RowEnd = size(digit,1);                       
bCheckEmpty = 0;

ce_count = 0;
for i=1:length(cc)
    if cc(i) == 0
        ce_count = ce_count+1;
        if ce_count > dist_empty
            bCheckEmpty = 1;
        end
    else
        ColumnStart = i-dist_empty;
        break;
    end
end

ce_count = 0;
for i=length(cc):-1:1
    if cc(i) == 0
        ce_count = ce_count+1;
        if ce_count > dist_empty
            bCheckEmpty = 1;
        end
    else
        ColumnEnd = i+dist_empty;
        break;
    end
end

ce_count = 0;
for i=1:length(rr)
    if rr(i) == 0
        ce_count = ce_count+1;
        if ce_count > dist_empty
            bCheckEmpty = 1;
        end
    else
        RowStart = i-dist_empty;
        break;
    end
end  

ce_count = 0;
for i=length(rr):-1:1
    if rr(i) == 0
        ce_count = ce_count+1;
        if ce_count > dist_empty
            bCheckEmpty = 1;
        end
    else
        RowEnd = i+dist_empty;
        break;
    end
end

if bCheckEmpty
    ColumnWidth = ColumnEnd-ColumnStart;
    RowHeight = RowEnd-RowStart;    
    digit = imcrop(digit,[ColumnStart RowStart ColumnWidth RowHeight]);
end

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
    
    pyth_y = max(5,round(.8*((ColumnEnd-ColumnStart)*thet/2)));
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
        
        set(line([ColumnStart, ColumnEnd], [RowStart, RowStart]), 'Color', box_color, 'Linewidth', 2);        
        set(line([ColumnStart, ColumnEnd], [RowEnd, RowEnd]), 'Color', box_color, 'Linewidth', 2);        
        set(line([ColumnStart, ColumnStart], [RowStart, RowEnd]), 'Color', box_color, 'Linewidth', 2);        
        set(line([ColumnEnd, ColumnEnd], [RowStart, RowEnd]), 'Color', box_color, 'Linewidth', 2);      
        
        center_row = out_pos(k,1);
        center_col = out_pos(k,2);
        
        set(line([center_col, center_col],[(center_row-7), (center_row+7)]),'color',center_color,'linewidth',1);
        set(line([(center_col-7), (center_col+7)],[center_row, center_row]),'color',center_color,'linewidth',1);
        drawnow;
%         disp(['Speed limit signboard detected -> score : ',num2str(out_val(k)),' aggregate : ',num2str(out_agg(k))]);        
    end  
else
%     disp(['No Circle Detected!']);
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