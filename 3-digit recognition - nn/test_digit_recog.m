function test_digit_recog

close all;
clear;

%-------------------------------%
nn_threshold = .5;
%-------------------------------%
load slr_data;
load slr_nn;
% load slr_nn2;
N = length(slr_sdata);

for i=1:N    
    numdigit = size(slr_sdata(i).num_ind,1);
    map = .5*ones(slr_height+20,(numdigit*slr_width+20));
    map(11:(10+slr_height),11:(numdigit*slr_width+10)) = slr_sdata(i).data;
    imshow(map); 
    str = ['image ',num2str(i)]; title(str);
    
    data = [];
    for j=1:numdigit
        temp = slr_sdata(i).data(:,slr_sdata(i).num_ind(j,1):slr_sdata(i).num_ind(j,2));
        data = [data temp(:)];
    end
    [output,str] = digit_recog(net,data,nn_threshold);    
    
    if output
        disp(['RESULT for image ',num2str(i),' : ',str,' km/h']);
    else
        disp(['RESULT for image ',num2str(i),' : ',str]);
    end
    pause(1);
end
disp(['completed...']);
