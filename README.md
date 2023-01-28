# 2022-PR-ERCAN

Using the code, please cite:
Wang J, Ma Z, Nie F, et al. Entropy Regularization for Unsupervised Clustering with Adaptive Neighbors[J]. Pattern Recognition, volume 125, 2022, 108517, doi: 10.1016/j.patcog.2021.108517.

https://www.sciencedirect.com/science/article/pii/S0031320321006932

The code explanation: 
The main function of the code: ERCAN_norm00k.m

If you have any questions, please connect zhenyu.ma@mail.nwpu.edu.cn

# Use of Main Function

[result,labelnew,r,t,ddd,S_sym] = ERCAN_norm00k(X,label,k,gamma,eta,rate,iter)

example: [result,~] = ERCAN_norm00k(X,label,10,1,1.05)

Input:
X: data matrix
label: ground truth (for clustering result)
k: the number of neighbors of similarity graph S
gamma: coefficient gamma (if nargin=4, gamma should be given by user)
eta: the multiplier of gamma (if nargin>4, gamma is computed by \eta multiples ||dij_x||_{F}/n, and gamma can be given any value)
rate: the change rate of Lambda, default 2
iter: Maximum of iteration, default 40
