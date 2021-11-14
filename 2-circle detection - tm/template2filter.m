function template2filter

%------------------------------------------------------%
img = imread('circle_template3.jpg');
[cir_height,cir_width,dummy] = size(img);
scale_ratio = 1.1;
fft_res = 250;
num_res = round(log10(fft_res/(cir_height+5))/log10(scale_ratio)+1);
filter_name = ['circle_filter'];
%------------------------------------------------------%

if isrgb(img)
    img = rgb2gray(img);
end
img = im2double(img);

new_img = img-mean(img(:));
new_img = imrotate(new_img,180);

for k=1:num_res
    sres = round(fft_res/(scale_ratio^(k-1)));
    temp = fft2(new_img,sres,sres);
    fft_cir(k,1) = {temp};
    clear temp;
end
 
str = ['save ',filter_name,' fft_cir cir_width cir_height scale_ratio;'];
eval([str]);
disp(['completed...']);