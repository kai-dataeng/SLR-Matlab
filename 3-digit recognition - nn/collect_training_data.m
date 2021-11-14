function collect_training_data

clear;
close all;

cur_dir = pwd;
cd(['digits\']);
file_img = dir('*.jpg');
N = length(file_img);
disp(['total images : ',num2str(N)]);


f1 = figure;
f2 = figure;
set(0,'Units','pixels'); 
scnsize = get(0,'ScreenSize');
bw = round(scnsize(3)/2);
ori_fsize = get(f2,'position');
bh = round((ori_fsize(4)/ori_fsize(3))*bw);
set(f2,'position',[-100 (scnsize(4)-bh-70) bw bh]);

%++++++++++++++++++++++++++++++++++++++++++++%
std_weight = 2;
min_obj_hscale = 1.15;
slr_sdata = [];
slr_width = 15;
slr_height = 25;
std_size = 100;
%++++++++++++++++++++++++++++++++++++++++++++%

for cnt = 1:N
    disp(['filename : ',file_img(cnt).name]);

    img = imread(file_img(cnt).name);
    img = imresize(img,[std_size std_size]);
    [row,col,dummy] = size(img);
    if dummy==3
        img = rgb2gray(img);
    end
    figure(f1),subplot(1,2,1),imshow(img);
    title(['original image']); drawnow;
    
    BW = im2bw(img,.5);
%     figure(f1),subplot(1,2,2),imshow(BW);
%     title(['binary images']); drawnow;

    BW_label = bwlabel(BW);
    num_obj = max(BW_label(:));
    temp_map = zeros(row,col);
    col_number = [1:col]; row_number = [1:row]';
    col_mat = repmat(col_number,row,1); row_mat = repmat(row_number,1,col);
    obj_ind = []; obj_area = []; digit_data = [];
    obj_count = 0;
    
    if num_obj>0
        for k=1:num_obj
            [r,c] = find(BW_label == k);
            obj_area = [obj_area; length(r)];
            obj_ind = [obj_ind; struct('row_ind',r,'col_ind',c)];
        end
        
        min_area = .3*max(obj_area);        
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
                
                ColumnStart = max(1, cMean - min(colStd, ceil(col/2)));
                ColumnEnd = min(col, cMean + min(colStd, ceil(col/2)));
                RowStart = max(1, rMean - min(rowStd, ceil(row/2)));
                RowEnd = min(row, rMean + min(rowStd, ceil(row/2)));                  
                
                figure(f1),subplot(1,2,2),imshow(temp_map); title(['detected digit']);                
                set(line([ColumnStart, ColumnEnd], [RowStart, RowStart]), 'Color', 'red', 'Linewidth', 2);        
                set(line([ColumnStart, ColumnEnd], [RowEnd, RowEnd]), 'Color', 'red', 'Linewidth', 2);        
                set(line([ColumnStart, ColumnStart], [RowStart, RowEnd]), 'Color', 'red', 'Linewidth', 2);        
                set(line([ColumnEnd, ColumnEnd], [RowStart, RowEnd]), 'Color', 'red', 'Linewidth', 2);              
                drawnow;                 
                
                ColumnWidth = max(1,(ColumnEnd-ColumnStart));
                RowHeight = max(1,(RowEnd-RowStart));
                obj_hscale = RowHeight/ColumnWidth;
                disp(['hscale digit: ',num2str(obj_hscale)]);
                
                if obj_hscale>min_obj_hscale
                    
                    obj_count = obj_count+1;
                    
                    if obj_hscale>=3
                        rStd = rStd*(std_weight+.1);
                        cStd = cStd*(std_weight+10);                        
                        ColumnStart = max(1, cMean - min(cStd, ceil(col/2)));
                        ColumnEnd = min(col, cMean + min(cStd, ceil(col/2)));
                        RowStart = max(1, rMean - min(rStd, ceil(row/2)));
                        RowEnd = min(row, rMean + min(rStd, ceil(row/2))); 
                        ColumnWidth = max(1,(ColumnEnd-ColumnStart));
                        RowHeight = max(1,(RowEnd-RowStart));
                        digit = imcrop(temp_map,[ColumnStart RowStart ColumnWidth RowHeight]);    
                    else                        
                        digit = imcrop(temp_map,[ColumnStart RowStart ColumnWidth RowHeight]);                    
                        digit = remove_empty_background(digit);
                    end
                    
                    figure(f2),imshow(digit);
                    title(['extracted digit']); drawnow;
                    
                    digit = imresize(digit,[slr_height slr_width]);
                    digit_data = [digit_data; struct('data',digit)];                    
                end                
                pause(.5);
                figure(f2),clf;                
                temp_map = temp_map*0;
            end
        end
        
        if obj_count>0            
            temp = []; temp2 = []; st = 1;
            for i=1:obj_count
                temp = [temp digit_data(i).data];
                temp2 = [temp2; st (slr_width+st-1)];
                st = i*slr_width+1;
            end                
            slr_sdata = [slr_sdata; struct('num_ind',temp2,'data',temp)];
        end            
    end    
    figure(f1),clf;
end
close all;

cd(['nondigits\']);
file_img = dir('*.jpg');
N = length(file_img);
disp(['read non-digit images']);
nr_sdata = [];
for cnt = 1:N
    disp(['filename : ',file_img(cnt).name]);
    img = imread(file_img(cnt).name);
    img = imresize(img,[slr_height slr_width]);
    [row,col,dummy] = size(img);
    if dummy==3
        img = rgb2gray(img);
    end
    img = im2double(img);
    nr_sdata = [nr_sdata; struct('data',img)]; 
end   
    
cd(cur_dir);
save slr_training_data slr_sdata nr_sdata slr_width slr_height;
disp(['completed...']);

%------------------------------------------------------------------------------------------%
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
    