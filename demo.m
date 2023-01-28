DataName = 'Coil20';
load([DataName,'_data.mat'])

k = 3;
gamma = 1;
eta = 10;
% —————— X: data matrix n*d
% —————— label: ground truth (for clustering result) n*1
% —————— k：neighbors of similarity matrix S

[result,~] = ERCAN(X,label,k,gamma,eta);