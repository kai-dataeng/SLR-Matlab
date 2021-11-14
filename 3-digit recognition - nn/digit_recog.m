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
