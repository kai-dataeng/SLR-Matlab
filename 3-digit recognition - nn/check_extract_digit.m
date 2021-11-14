function check_extract_digit

close all;
clear;

load slr_data;
% load slr_training_data;
N = length(slr_sdata);

for i=1:N    
    numdigit = size(slr_sdata(i).num_ind,1);
    map = .5*ones(slr_height+20,(numdigit*slr_width+20));
    map(11:(10+slr_height),11:(numdigit*slr_width+10)) = slr_sdata(i).data;
    imshow(map); 
    str = ['image ',num2str(i)];
    title(str);
    pause(2);
end
disp(['completed...']);
