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
S2 = num_class;
N = max_N;
data_range = [zeros((slr_width*slr_height),1) ones((slr_width*slr_height),1)];
net = newff(data_range,[S1 S2],{'logsig' 'logsig'},'trainscg');
net = init(net);

net.performFcn = 'sse';        % Sum-Squared Error performance function
net.trainParam.goal = 5e-6;     % Sum-squared error goal.
net.trainParam.show = 10;      % Frequency of progress displays (in epochs).
net.trainParam.epochs = 5e3;  % Maximum number of epochs to train.
net.trainParam.lr = 0.1;
net.trainParam.mc = 0.95;      % Momentum constant.

P = []; T = [];
target = eye(S2);

ind_nd = [];
for j=1:N
    for i=1:num_class
        temp = slr_sdata(i).data(:,slr_sdata(i).num_ind(j,1):slr_sdata(i).num_ind(j,2));
        P = [P temp(:)];
    end
    T = [T target];
end
net = train(net,P,T);

disp(['NN training completed..']);
save slr_nn net slr_height slr_width;
validate_nn(net,P,T);