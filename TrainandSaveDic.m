function [D,energy, Etime] = TrainandSaveDic(Z,K,T,patchnum,lamda)
%%                  input
% Z                is the patch set that used to trained the dictionary, 
% K                is the number of atom in the dictionary, 
% T                is the iterations we trained the dictionary
% patchnum         is how many patches used in one iteration 
% lamda            is the parameter that in the SUNSAL algorithm

%%                 output
% D                the trained dictionary

param.K=K;  % learns a dictionary with 
param.T=T;  % let us see what happens after 500 iterations.
param.patchnum = patchnum;
param.lamda = lamda;
param.err = 1e-3;

%%%%%%%%%% FIRST EXPERIMENT %%%%%%%%%%%

% K = 256;
% T = 500;
% p = 64;
% lamda = 1.8;


% tic;
% [D] = DicLearning(X,K,T,p,lamda);
% toc
patchsize = sqrt(size(Z,1));

[D, energy, Etime] = DicLearningM(Z,param);
savefile = sprintf('dic%dX%d.mat',patchsize,patchsize);

save(savefile,'D');
% image_esti = (image_esti-min(image_esti(:)))/(max(image_esti(:))-min(image_esti(:)));