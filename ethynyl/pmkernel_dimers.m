clearvars;
close all;
clc;

LIBSVMPATH = '~/Documents/ML/libsvm-3.23/matlab';

addpath('svm');
addpath(LIBSVMPATH);

datapath = './CarboDimersVdwP2.mat';
load(datapath);


% change cell to struct
%new = cell2struct(AM(:,1)',{'am'})


% run kernel thing
K = pmkernel_unlabeled(new, 4, 6); %weights for kernel

Knew_scaled = rescale(K);
% 0 in diagonals
Knew_scaled(logical(eye(size(Knew_scaled)))) = 0;


csvwrite('Kernel_sparse_VdwP2.csv',Knew_scaled)