function slr_initialize

load('D:\WORK\work - Matlab II\speed limit recognition\2-circle detection - tm\circle_filter');
min_res = 600;
h_end = 6;
h_sta = 4;
h_bin = 9;
v_end = 5;
v_sta = 1;
v_bin = 6;
resolution = 500;
theta = -pi/9;
box_color = ['green'];
center_color = ['blue'];
num_scales = size(fft_cir,1);
fft_res = size(fft_cir{1},1);

% cir_threshold = 0.7;
% max_corr = 150;  
cir_threshold = 0.4;
max_corr = 80;  

da_size = 100;
uw_area = round(da_size*da_size*1e-2);
temp_map = zeros(da_size,da_size);
col_number = [1:da_size]; row_number = [1:da_size]';
col_mat = repmat(col_number,da_size,1); row_mat = repmat(row_number,1,da_size);
std_weight = 2;
min_obj_hscale = 1.1;
nn_threshold = .5;

save slr_init_data