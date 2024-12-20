function result = main_fuction(anchor_rate,order,k,beta,rng_seed,X,Y)
%% load data
% 1-Glass; 2-Yeast; 3-Wine; 4-German; 5-Dermatology
% 6-JAFFE; 7-COIL-20; 8-MSRA-25
[n, dim] = size(X);
dimY = length(Y);
c = length(unique(Y));

if n ~=dimY
    X = X';
    dim = n;
    n = dimY;
end 
%% data pre-process
X = X./max(X,[],2);% n x dim

%% hyperparameter setting
anchor_rate = anchor_rate;
maxIter = 50;
opts.anchorSelectStyle = 3;
order = order;
beta = beta;
k = k;
rng(rng_seed);

Isconverg = 0;
eta = 1.1;
mu = 10e-5;
rho = 10e-5;
max_mu = 10e12;
max_rho = 10e12;


m = fix(n*anchor_rate);
%%
% disp('----------Anchor Selection----------');
if opts.anchorSelectStyle == 1 % direct sample
    [~,ind,~] = graphgen_anchor(X,m);
    centers = X(ind, :);% m x dim
elseif opts.anchorSelectStyle == 2 % rand sample
    vec = randperm(n);
    ind = vec(1:m);
    centers = X(ind, :);
elseif opts.anchorSelectStyle == 3 % KNP
    [~, ~, ~, ~, dis] = litekmeans(X,m);
    [~,ind] = min(dis,[],1);
    ind = sort(ind,'ascend');
    centers = X(ind, :);
elseif opts.anchorSelectStyle == 4 % kmeans sample
    [~, centers, ~, ~, ~] = litekmeans(X, m);
end
%%
% disp('----------1st order 2P Graphs Inilization----------');
D = L2_distance_1(X', centers'); % n x m
[~, idx] = sort(D, 2);
B1 = zeros(n,m);
for ii = 1:n
    id = idx(ii,1:k+1);
    di = D(ii, id);
    B1(ii,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
%%
% disp('----------Generate high order 2P Graphs----------');
B = cell(order,1); % cell;order x 1 -> n x m
B{1} = B1./max(max(B1,[],2));
[U,sigma,Vt] = svd(B1);
for d = 2:order
    temp = U*sigma.^(2*d-1)*Vt';
    temp(temp<eps)=0;
    temp = temp./max(max(temp,[],2));
    %temp = temp./sum(temp,2);
    B{d,1} = temp;
end

%% initial variable
B_hat = time2frequency(B);
[initial_y, ~, ~, ~, dis] = litekmeans(X,c);
F_init_ = n2nc(initial_y ,c);
% indices = randperm(length(initial_y));
% initial_y = initial_y(indices);
% F_init_ = n2nc(initial_y ,c);
G_init_ = n2nc(initial_y(ind) ,c);
for o = 1:d 
    F_init{o} = F_init_;
    G_init{o} = G_init_;
end 
F_hat = time2frequency(F_init);
G_hat = time2frequency(G_init);

for o = 1:d
    Y1_r{o} = zeros(n, c);   % about F low rank         
    Q_r{o} = zeros(n, c);    % about F low rank          
    Y2{o} = zeros(n, c);   % about F > 0
    J{o} = zeros(n, c);    % about F > 0
end
iter = 1;
% % Y1
% Y1_r_tensor = cat(3,Y1_r{:});
% Y1_tensor = permute(Y1_r_tensor,[1,3,2]);
% for o = 1:d
%     Y1{o} = Y1_tensor(:,:,o);
% end
% % Q
% Q_r_tensor = cat(3,Q_r{:});
% Q_tensor = permute(Q_r_tensor,[1,3,2]);
% for o = 1:d
%     Q{o} = Q_tensor(:,:,o);
% end

tic;
while(Isconverg == 0) 
    fprintf('This is iter: %d \n', iter);
    %% 1. update G_hat{v}
    for o = 1:d
        G_hat{o} = B_hat{o}' * F_hat{o};
    end
    G = frequency2time(G_hat);
    
    %% 2. update F_hat{v}
    Q_r_hat = time2frequency(Q_r);
    Y1_r_hat = time2frequency(Y1_r);
    J_hat = time2frequency(J);
    Y2_hat = time2frequency(Y2);
    
    for o = 1:d
        W_hat{o} = B_hat{o} * G_hat{o} + mu * Q_r_hat{o} + Y1_r_hat{o} + rho * J_hat{o} + Y2_hat{o};
        [nn{o}, ~, vv{o}] = svd(W_hat{o}, 'econ');
        F_hat{o} = nn{o} * vv{o}';
    end
    F = frequency2time(F_hat);
    
    %% 3. update Q_r{v}
    F_tensor = cat(3,F{:});
    F_r_tensor = permute(F_tensor,[1,3,2]);
    Y1_r_tensor = cat(3,Y1_r{:});
    Y1_tensor = permute(Y1_r_tensor,[1,3,2]);
    M_tensor = F_r_tensor - Y1_tensor./mu;
    [Q_tensor,tnn,trank] = solve_A_tensor(M_tensor,beta/mu);
    Q_r_tensor = permute(Q_tensor,[1,3,2]);
    for o = 1:d
        Q_r{o} = Q_r_tensor(:,:,o);
    end
    Q_r_hat = time2frequency(Q_r);
    
    
    %% 4. update  J{v}
    for o = 1:d
        J{o} = F{o} + Y2{o} ./ rho;
        J{o}(J{o}<0) = 0;
    end
    J_hat = time2frequency(J);
    
    %% update Y1 Y2 mu rho
    
    for o = 1:d
        Y1_r{o} = Y1_r{o} + mu * (F{o} - Q_r{o});
        Y2{o} = Y2{o} + rho * (F{o} - J{o});
    end
    Y1_r_hat = time2frequency(Y1_r);
    Y2_hat = time2frequency(Y2);
    mu = min(eta * mu, max_mu);
    rho = min(eta * rho, max_rho); 
    
    %%
    if iter == maxIter
        Isconverg = 1;
    end
    %% 
    for o = 1:d
        T{o} = B_hat{o} - F_hat{o} * G_hat{o}';
        nm(o) = norm(T{o},'fro');
    end
    obj_fun(iter) = sum(nm) + tnn;
    if iter > 1 && abs(obj_fun(iter)-obj_fun(iter-1)) < 0.001
        Isconverg = 1;
    end
    
    
    iter = iter + 1;     
end
time_record = toc;
F_sum = F{1};
for o = 2:d
    F_sum = F_sum + F{o};
end
F_final = F_sum / d;
[~, Y_pre] = max(F_final, [], 2); 
my_result = ClusteringMeasure1(Y, Y_pre);
fprintf('acc = %f !!!!! \n',my_result(1));
result = [my_result time_record];
end