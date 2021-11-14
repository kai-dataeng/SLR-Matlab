function validate_nn(net,data,label)

Y = sim(net,data);

Ptarget = max(label(:));
Ntarget = min(label(:));
threshold = (Ptarget+Ntarget)/2;

N = size(Y,2);

TAR = 0; FAR = 0;
for i=1:N
    Tind = find(label(:,i) == 1);
    [max_score,Pind] = max(Y(:,i));
    if max_score >= threshold
        if Tind == Pind
            TAR = TAR+1;
        else
            FAR = FAR+1;
        end
    end
end
        
FAR = FAR/N;
TAR = TAR/N;

disp(['-------------------------------------------------------']);
disp(['False Accepted Rate (FAR) : ',num2str(FAR*100),'%']);
disp(['True Accepted Rate (TAR) : ',num2str(TAR*100),'%']);
disp(['-------------------------------------------------------']);
