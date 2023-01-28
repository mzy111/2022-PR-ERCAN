function [S_init,gamma,X_exp,sum_i_exp] = init_S(X,gamma,eta)
% X input data, n*d
% S weight matrix
% min sij||x_i-x_j||^2-gamma(sij*logsij), s.t.S'1=1, 0=<s_ij<=1  (gamma>=0��gamma��ֵ��Ӧ��ͬ���ݼ�Ҫȡ��һ������ֵ��ͨ���ܴ�)
% use Lagrange multipliers to initialize Sij
% S_initΪ��ʼ����Ȩ�ؾ���S
% gammaΪʹX_e���󲻳���inf�Ľ���ֵ�������ţ�
% E_initΪ��Ӧ��ʼS�ĳ�ʼ�أ�Ŀ�꺯����
% pָ��gamma�����ѡȡ��������F����ȷ��


[n,~] = size(X);                    % ������Ŀ
X_exp = zeros(n,n);                 % �ݾ����м����
S_init = zeros(n,n);                % ����S
dij_x = L2_distance_1(X',X');       % ������������������
if nargin == 3
    gamma = norm(dij_x,'fro')/n;        % ����F��������Ԥ��ʼ��gamma
    gamma = eta*gamma;                % ��������
elseif nargin == 1
    error('The parameter gamma should be input before initialization.')
end


%% ��ʼ������
sum_i_exp = zeros(n,1);  % ��ѭ����ʱ�临�Ӷ�Ϊn^{2}
for i = 1:n
    for j = 1:n
        X_exp(i,j) = exp(-(1/gamma)*dij_x(i,j));  % ��gamma��˲�����   �����i�е����ƶ���
    end
    sum_i_exp(i) = sum(X_exp(i,:));               % ������еĺ�            
    for j = 1:n
        S_init(i,j) = X_exp(i,j)/(sum_i_exp(i));  % �����i�е����ƶ�
    end
end