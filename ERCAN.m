function [result,labelnew,r,t,ddd,S_sym] = ERCAN(X,label,k,gamma,eta,rate,iter)

% XΪԭʼ����n*d
% labelΪ��ʵ��ǩ����n*1��������֤ERCAN�������ľ���ָ��
% �÷�������||s_i||0 = k�����ƶȾ������ǿϡ��Լ����kΪͼ������
% �������Ϊ5����Ĭ��gamma��F����ȷ����sigma������������
% �������Ϊ4����Ĭ��gamma��Ϊ������sigma��ֵΪ0
% �������Ϊ3���������ѱ���δ����gammaֵ
% iterΪ����������
% rateΪlambda�ĵ������ʣ�Ĭ��Ϊ2

% result: ����ACC,NMI,Purity��Running Time��4ά������
% labelnew: ERCAN�㷨������n*1
% r: Rank of Laplacian matrix
% t: Running Time
% ddd: ��ǰѭ������
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
% self_similar 0��ÿ��������������������ƶȱ��޶�Ϊ0����0��ÿ����������������������ƶ�
% k            0����Լ��������||s_i||0 = k,��self_similar = 1,��||s_i||0 = k+1



tic
pause(1)
%% Parameter Setting
[n,~] = size(X);
c = length(unique(label));   % the number of class
dij_x = L2_distance_1(X',X');% ƽ���������E n*n ʱ�临�Ӷ�Ϊn^{2}*d





%% Step1����ʼ��S0����F0����
if eta == 0
    [S,gamma,~,~] = init_S(X,gamma);        % S�����gamma�ĳ�ʼ����S0��gamma0
else
    [S,gamma,~,~] = init_S(X,gamma,eta);  % ʱ�临�Ӷ�Ϊn^{2}*d + n^{2}
end
% S = rand(n,n);
% sumS = sum(S,2);
% S = S./sumS;
% gamma = gamma*sigma;
Lambda = gamma;                                     % ��ʼ��lambda��ֵ��������˹������Լ�����򻯲���
S_init = (S + S')/2;                                % �Գƻ�
D_0 = diag(sum(S_init,2));                          % ��Ⱦ���
L_0 = D_0 - S_init;                                 % ��������˹����
[F,~,~] = eig1(L_0, c, 0);                          % ѡǰc�е����������Լ���С��c������ֵ ʱ�临�Ӷ�Ϊn^{3}

r0 = rank(L_0);
if r0 < n-c                                  % �б���ͨ�����Ĵ�С,��ͨ����Ϊ1˵��ȫ����ͨ����ͨ��������˵��ƾ�ճ����˲����ڵ������Ҫ����
    error('The original graph has more than %d connected component', c);
end




%% Step2���̶�F������S  ʱ�临�Ӷ�Ϊn^{2}
sum_i_original = zeros(n,1);
% dia = zeros(n,1);
for ddd = 1:iter  
    dij_f = L2_distance_1(F',F');                                     % ��dij_f
    dij = zeros(n,n);                                                 % dij��ʼ
    D_exp_update = zeros(n,n); 
    % ���ݾ����ʼ
    for i = 1:n
        for j = 1:n
            dij(i,j) = dij_x(i,j) + Lambda*dij_f(i,j);                % ��dij
            D_exp_update(i,j) = exp(-(1/gamma)*dij(i,j));             % ���ݾ���
        end
%         dia(i) = D_exp_update(i,i);
        sum_i_original(i) = sum(D_exp_update(i,:));
        [~,idx] = sort(D_exp_update(i,:),'descend');                  % �Ӵ�С��������
        id = idx(k+2:n);                                              % ѡ������k+1��n�ľ���Ԫ��
        D_exp_update(i,id) = 0;                                       % �ñȽϴ�����Ԫ��ȫ��Ϊ0
        if self_similar == 0
            D_exp_update(i,i) = 0;
        end
        sum_i_update = sum(D_exp_update(i,:));                        % �õ�i��k��Ԫ�ص��к�
        for j = 1:n
            S(i,j) = D_exp_update(i,j)/(sum_i_update);                % ��ÿ��Ԫ�ض��к͵�ռ�ȣ�sij�Ĺ�ʽ�Ƶ����
        end
    end
    

%% ���������̶�S������F ʱ�临�Ӷ�Ϊn^{3}

    S_sym = (S+S')/2;               % �Գƻ�                       
    D_s = diag(sum(S_sym,2));       % �Ⱦ���
    L_s = D_s - S_sym;              % ������˹����
    r = rank(L_s);                  % ��¼������˹�������
    F_old = F;                      % ����F����
    [F,~,~]=eig1(L_s, c, 0);        % ѡǰc�е����������Լ���С��c������ֵ,ev�ǰ���n������ֵ��������
    ddd
    %evs(:,ddd+1) = ev;             % evs�ĵ�i�ж�Ӧ��i�ε���������ֵ������


%% ��������lamda
%     fn1 = sum(ev(1:c));
%     fn2 = sum(ev(1:c+1)); 
    %if fn1 > 0.00000000001     % ���fn1�Ƚϴ�˵��connected componentsС��c����Ҫ����lambda
    [cluster, ~] = graphconncomp(sparse(S_sym));
    if cluster < c              % cc̫�٣���Ҫ��ǿ��Լ����ǿ��
        Lambda = Lambda*rate;   % ����cc
    %elseif fn2 < 0.00000000001 % ���fn1��С��˵��connected components����c����Ҫ����lambda
    elseif cluster > c          % cc̫����Ҫ��С��Լ����ǿ�ȣ���������F����S���ص���75��
        Lambda = Lambda/rate;  F = F_old;
    else
        break;
    end
end
[~, labelnew] = graphconncomp(sparse(S_sym)); 
labelnew = labelnew';           % ������˹������ж��Ƕ�Ӧ������������ƶȾ��������Ҫ�б�ľ����ǶԳƻ�֮���S����
t = toc;
result = zeros(1,4);
result(1:3) = ClusteringMeasure(label,labelnew);
result(4) = t;