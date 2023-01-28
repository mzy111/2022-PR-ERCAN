function [result,labelnew,r,t,ddd,S_sym] = ERCAN(X,label,k,gamma,eta,rate,iter)

% X为原始数据n*d
% label为真实标签向量n*1，用于验证ERCAN聚类结果的聚类指标
% 该方法利用||s_i||0 = k对相似度矩阵进行强稀疏约束，k为图近邻数
% 如果输入为5个，默认gamma由F范数确定，sigma用来调整倍率
% 如果输入为4个，默认gamma人为给定，sigma的值为0
% 如果输入为3个，会提醒报错，未给定gamma值
% iter为最大迭代次数
% rate为lambda的调整倍率，默认为2

% result: 包含ACC,NMI,Purity和Running Time的4维行向量
% labelnew: ERCAN算法聚类结果n*1
% r: Rank of Laplacian matrix
% t: Running Time
% ddd: 当前循环次数
% S_sym: Symmetric Similarity matrix

self_similar = 0;
if nargin < 7
    iter = 40;
end
if nargin < 6
    rate = 2;
end
if nargin < 5
    eta = 0;
end
if nargin == 3
    error('The parameter gamma should be input.')
end
% self_similar 0：每个样本自身与自身的相似度被限定为0；非0：每个样本自身与自身存在相似度
% k            0范数约束参数：||s_i||0 = k,若self_similar = 1,则||s_i||0 = k+1



tic
pause(1)
%% Parameter Setting
[n,~] = size(X);
c = length(unique(label));   % the number of class
dij_x = L2_distance_1(X',X');% 平方距离矩阵E n*n 时间复杂度为n^{2}*d





%% Step1：初始化S0，求F0矩阵
if eta == 0
    [S,gamma,~,~] = init_S(X,gamma);        % S矩阵和gamma的初始化：S0和gamma0
else
    [S,gamma,~,~] = init_S(X,gamma,eta);  % 时间复杂度为n^{2}*d + n^{2}
end
% S = rand(n,n);
% sumS = sum(S,2);
% S = S./sumS;
% gamma = gamma*sigma;
Lambda = gamma;                                     % 初始化lambda的值：拉普拉斯矩阵秩约束正则化参数
S_init = (S + S')/2;                                % 对称化
D_0 = diag(sum(S_init,2));                          % 求度矩阵
L_0 = D_0 - S_init;                                 % 求拉普拉斯矩阵
[F,~,~] = eig1(L_0, c, 0);                          % 选前c列的特征向量以及最小的c个特征值 时间复杂度为n^{3}

r0 = rank(L_0);
if r0 < n-c                                  % 判别连通分量的大小,连通分量为1说明全体连通，连通分量过多说明凭空出现了不存在的类别，需要报错
    error('The original graph has more than %d connected component', c);
end




%% Step2：固定F，更新S  时间复杂度为n^{2}
sum_i_original = zeros(n,1);
% dia = zeros(n,1);
for ddd = 1:iter  
    dij_f = L2_distance_1(F',F');                                     % 求dij_f
    dij = zeros(n,n);                                                 % dij初始
    D_exp_update = zeros(n,n); 
    % 求幂矩阵初始
    for i = 1:n
        for j = 1:n
            dij(i,j) = dij_x(i,j) + Lambda*dij_f(i,j);                % 求dij
            D_exp_update(i,j) = exp(-(1/gamma)*dij(i,j));             % 求幂矩阵
        end
%         dia(i) = D_exp_update(i,i);
        sum_i_original(i) = sum(D_exp_update(i,:));
        [~,idx] = sort(D_exp_update(i,:),'descend');                  % 从大到小距离排序
        id = idx(k+2:n);                                              % 选出其余k+1到n的距离元素
        D_exp_update(i,id) = 0;                                       % 让比较大距离的元素全部为0
        if self_similar == 0
            D_exp_update(i,i) = 0;
        end
        sum_i_update = sum(D_exp_update(i,:));                        % 该第i行k个元素的行和
        for j = 1:n
            S(i,j) = D_exp_update(i,j)/(sum_i_update);                % 求每个元素对行和的占比：sij的公式推导结果
        end
    end
    

%% 步骤三：固定S，更新F 时间复杂度为n^{3}

    S_sym = (S+S')/2;               % 对称化                       
    D_s = diag(sum(S_sym,2));       % 度矩阵
    L_s = D_s - S_sym;              % 拉普拉斯矩阵
    r = rank(L_s);                  % 记录拉普拉斯矩阵的秩
    F_old = F;                      % 储存F矩阵
    [F,~,~]=eig1(L_s, c, 0);        % 选前c列的特征向量以及最小的c个特征值,ev是包含n个特征值的列向量
    ddd
    %evs(:,ddd+1) = ev;             % evs的第i列对应第i次迭代的特征值列向量


%% 调整参数lamda
%     fn1 = sum(ev(1:c));
%     fn2 = sum(ev(1:c+1)); 
    %if fn1 > 0.00000000001     % 如果fn1比较大，说明connected components小于c，需要增大lambda
    [cluster, ~] = graphconncomp(sparse(S_sym));
    if cluster < c              % cc太少，需要加强秩约束的强度
        Lambda = Lambda*rate;   % 增大cc
    %elseif fn2 < 0.00000000001 % 如果fn1很小，说明connected components大于c，需要减少lambda
    elseif cluster > c          % cc太大，需要减小秩约束的强度，并重新用F更新S，回到第75行
        Lambda = Lambda/rate;  F = F_old;
    else
        break;
    end
end
[~, labelnew] = graphconncomp(sparse(S_sym)); 
labelnew = labelnew';           % 拉普拉斯矩阵的判定是对应于组成它的相似度矩阵，因此需要判别的矩阵是对称化之后的S矩阵
t = toc;
result = zeros(1,4);
result(1:3) = ClusteringMeasure(label,labelnew);
result(4) = t;