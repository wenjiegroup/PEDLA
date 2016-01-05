function [ stackedAEOptTheta,netconfig,stackedAETheta ] = unsupervisedlearnig(train_x,sizes, flag, opts )
%DEEPLEARNIG Summary of this function goes here
%   Detailed explanation goes here
%
%train_x: m*n matrix. m, number of sample;n ,feature dimension of input
%train_y: m*1 vector. the value must be integer,and begin from 1 to the number of class .for
%           example 1,2,3,4,... 
%test_x,test_y:same as training set
%sizes: unit number of hidden layer. must be a vector. for example [200 200]
%opts: used to set numepochs, batchsize, momentum and alpha . you can use
%           the default

if nargin<2
    sizes=[200 200];
end
if nargin<3
    flag=0;
end
if ~exist('opts', 'var')
    opts = struct;
end
if ~isfield(opts, 'numepochs')   
    opts.numepochs = 1;
end
if ~isfield(opts, 'batchsize')   
    opts.batchsize = 100;
end
if ~isfield(opts, 'momentum')   
    opts.momentum = 0;
end
if ~isfield(opts, 'alpha')   
    opts.alpha = 1;
end

addpath('./minFunc/');
lambda = 3e-3;       % weight decay parameter  

rng(0);
ticid=tic;

%% unsupervised feature learning
dbn = dbnsetup(sizes, train_x, opts);
dbn = dbntrain(dbn, train_x, opts);


%% stack dbn
stack = cell(numel(sizes),1);
for i=1:numel(sizes)
    stack{i}.w=dbn.rbm{i}.W;
    stack{i}.b=dbn.rbm{i}.c;
end
for i=1:numel(sizes)
    stack{i+numel(sizes)}.w=dbn.rbm{numel(sizes)+1-i}.W';
    stack{i+numel(sizes)}.b=dbn.rbm{numel(sizes)+1-i}.b;
end
[stackparams, netconfig] = stack2params(stack);
stackedAETheta = [ stackparams ];

%% Finetune the whole model
%opts.numepochs = 50;
stackedAEOptTheta=unsupervisedfinetuning(stackedAETheta, train_x,...
                 netconfig, lambda, opts,flag );

T=toc(ticid);
fprintf('The total time: %.4f\n',T);

end

