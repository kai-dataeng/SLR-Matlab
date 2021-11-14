function slr_training

close all;
clear;

load slr_training_data;
num_class = length(slr_sdata);

max_N = 0;
for i=1:num_class
    cur_N = length(slr_sdata(i).num_ind(:,1));
    if cur_N > max_N
        max_N = cur_N;
    end
end
    
for i=1:num_class
    cur_N = length(slr_sdata(i).num_ind(:,1));
    if cur_N<max_N
        dd = ceil(max_N/cur_N);
        ind = [];
        count = 0;
        for j=1:dd
            permN = randperm(cur_N);
            for k=1:cur_N
                if count>(max_N-1)
                    break;
                end
                ind = [ind; slr_sdata(i).num_ind(permN(k),:)];
                count = count+1;
            end
        end
        slr_sdata(i).num_ind = ind;
    end
end

S1 = 15;
S2 = num_class+1;
N = max_N;
data_range = [zeros((slr_width*slr_height),1) ones((slr_width*slr_height),1)];
net = newff(data_range,[S1 S2],{'logsig' 'logsig'},'trainscg');
net = init(net);

net.performFcn = 'sse';        %Sum-Squared Error performance function
net.trainParam.goal = 5e-6;     %Sum-squared error goal.
net.trainParam.show = 10;      %Frequency of progress displays (in epochs).
net.trainParam.epochs = 5e3;  %Maximum number of epochs to train.
net.trainParam.lr = 0.1;
net.trainParam.mc = 0.95;      %Momentum constant.

P = []; T = [];
target = eye(S2);

ind_nd = [];
for j=1:N
    for i=1:num_class
        temp = slr_sdata(i).data(:,slr_sdata(i).num_ind(j,1):slr_sdata(i).num_ind(j,2));
        P = [P temp(:)];
    end
    P = [P zeros((slr_width*slr_height),1)];    
    ind_nd = [ind_nd ((num_class+1)*j)];
    T = [T target];
end

net = train(net,P,T);
disp(['Initial training NN completed..']);

%--------------------------------------%
noise_level = [.1:.05:.5];
Nnl = length(noise_level);
degree = [1:360];
Ndg = length(degree);
%--------------------------------------%

net.trainParam.goal = 1e-3;     % Sum-squared error goal.
net.trainParam.epochs = 1e3;  % Maximum number of epochs 
Pnew = P;

num_nd = length(nr_sdata);
for k=1:num_nd    
    nd_img = nr_sdata(k).data;
    nd_data = [];    
    for i=1:N
        rand_syth = randperm(8);        
        switch rand_syth(1)
            case 1
                n = randperm(Nnl);
                nd_temp = nd_img + randn(slr_height,slr_width)*noise_level(n(1));
            case 2
                nd_temp = fliplr(nd_img);
            case 3
                nd_temp = flipud(nd_img);            
            case 4
                n = randperm(Ndg);                
                nd_temp = imrotate(nd_img,degree(n(1)),'nearest','crop');
            case 5
                n = randperm(Nnl);
                nd_temp = nd_img + randn(slr_height,slr_width)*noise_level(n(1));
                nd_temp = fliplr(nd_temp);
            case 6
                n = randperm(Nnl);
                nd_temp = nd_img + randn(slr_height,slr_width)*noise_level(n(1));
                nd_temp = flipud(nd_temp);
            case 7
                n = randperm(Nnl);
                nd_temp = nd_img + randn(slr_height,slr_width)*noise_level(n(1));               
                n = randperm(Ndg);                
                nd_temp = imrotate(nd_temp,degree(n(1)),'nearest','crop');
            otherwise
                nd_temp = nd_img;
        end
        nd_temp = im2bw(nd_temp,.5);
        nd_data = [nd_data nd_temp(:)];         
    end
    Pnew(:,ind_nd) = nd_data;
    net = train(net,Pnew,T);
    disp(['training NN : ',num2str(k),' pass..']);
end

net.trainParam.goal = 1e-5;     % Sum-squared error goal.
net.trainParam.epochs = 5e3;  % Maximum number of epochs to train.
net = train(net,P,T);

disp(['NN training completed..']);
save slr_nn2 net slr_height slr_width;
validate_nn(net,P,T);