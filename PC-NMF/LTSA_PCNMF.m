function [output, W, W_cluster, H]=LTSA_PCNMF(data, basis_exchange_clusters, segment_width, basis_num, iter_num)
% This program aims to separate different sound sources with different levels of periodicity by using the PC-NMF.
% 
% Input parameters
% data: a non-negative spectrogram for source separation
% basis_exchange_clusters: number of source you want to separate, it can also be a range if the number of source is unknown
% segment_width: number of time frames you want to cascade the input spectrogram
% basis_num: number of basis in NMF
% iter_num: number of iterations in NMF
%
% output: spectrogram of each sound source separated by using PC-NMF
% W: basis matrix provided by NMF
% W_cluster: numerical label of each sound source for separating the basis matrix
% H: encoding matrix provided by NMF
%
% Written by Tzu-Hao Harry Lin (Academia Sinica, Taiwan), 2017/5/11
% Corresponding email: schonkopf@gmail.com

%% Parameters
if nargin<2
   basis_exchange_clusters=2; % number of clusters to unsupervised separate groups with different periodicity 
end
if nargin<3
    segment_width=1; % number of sequential vectors
end
if nargin<4
    basis_num=60; % number of basis for NMF
	iter_num=200; % number of iterations for updating the basis and encoding matrix during NMF
end

sparseness_W=[]; % sparsness parameter for updating the basis matrix during NMF (1 means very sparse, 0 means toward the gaussian distribution)
sparseness_H=[]; % sparsness parameter for updating the encoding matrix

%% Main Procedures
f_dim=size(data,1); 
data=sequential_matrix(data, segment_width);

% pretrain the NMF for 10 times and find the best initials
pretrain_iter_num=10;
for m=1:10
   [W(:,:,m),H(:,:,m)] = nmfsc(data,basis_num,sparseness_W,sparseness_H,pretrain_iter_num,0);
   V=W(:,:,m)*H(:,:,m);V=(V-min(min(V)))./max(max(V));
   rmse(m)=mean(reshape((data-V).^2,[],1))^0.5;
end
[~,c]=min(rmse);
[W,H] = nmfsc(data,basis_num,sparseness_W,sparseness_H,iter_num-pretrain_iter_num,0,W(:,:,c),H(:,:,c));

output=[]; clear V;
% Cluster different basis matrices with different periodicity by sparse NMF
[W, H, W_cluster,period_feature]=basis_exchange(W, H, basis_exchange_clusters, segment_width);

% Use traditional NMF enhancement procedure to reconstruct different sound sources
for m=1:max(W_cluster)
    output(:,:,m)=matrix_mean(data.*(W(:,W_cluster==m)*H(W_cluster==m,:)./(W*H)),segment_width,f_dim);
end


