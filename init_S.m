function [S_init,gamma,X_exp,sum_i_exp] = init_S(X,gamma,eta)
% X input data, n*d
% S weight matrix
% min sij||x_i-x_j||^2-gamma(sij*logsij), s.t.S'1=1, 0=<s_ij<=1  (gamma>=0且gamma的值对应不同数据集要取不一样的数值，通常很大)
% use Lagrange multipliers to initialize Sij
% S_init为初始化的权重矩阵S
% gamma为使X_e矩阵不出现inf的较优值（非最优）
% E_init为对应初始S的初始熵（目标函数）
% p指代gamma是随机选取还是利用F范数确定


[n,~] = size(X);                    % 样本数目
X_exp = zeros(n,n);                 % 幂矩阵（中间矩阵）
S_init = zeros(n,n);                % 矩阵S
dij_x = L2_distance_1(X',X');       % 基于样本构造距离矩阵
if nargin == 3
    gamma = norm(dij_x,'fro')/n;        % 利用F范数进行预初始化gamma
    gamma = eta*gamma;                % 调整倍率
elseif nargin == 1
    error('The parameter gamma should be input before initialization.')
end


%% 初始化过程
sum_i_exp = zeros(n,1);  % 该循环的时间复杂度为n^{2}
for i = 1:n
    for j = 1:n
        X_exp(i,j) = exp(-(1/gamma)*dij_x(i,j));  % 与gamma相乘并求幂   计算第i行的相似度行
    end
    sum_i_exp(i) = sum(X_exp(i,:));               % 计算该行的和            
    for j = 1:n
        S_init(i,j) = X_exp(i,j)/(sum_i_exp(i));  % 计算第i行的相似度
    end
end