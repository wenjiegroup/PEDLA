function [ acc,stackedAEOptTheta,output1,output2,output3, T ] = supervisedlearnig(train_x,train_y,test_x,test_y, sizes, flag, opts )
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

if nargin<5
    sizes=[200 200];
end
if nargin<6
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

%addpath('./minFunc/');
softmaxlambda=1e-4; % Weight decay parameter for softmax

% rng(0);
ticid=tic;

%% unsupervised feature learning
dbn = dbnsetup(sizes, train_x, opts);
dbn = dbntrain(dbn, train_x, opts);

%% input of softmax
trainfeatue=train_x;
for i=1:numel(sizes)
    trainfeatue=rbmup(dbn.rbm{i},trainfeatue);
end
inputSize=size(trainfeatue,2);
numClasses=length(unique(train_y));
softmaxModel=softmaxtrain(inputSize,numClasses , trainfeatue, train_y, opts, softmaxlambda,flag );
%softmaxoutput=softmaxtup( softmaxModel.optTheta , trainfeatue');
saeSoftmaxOptTheta = softmaxModel.optTheta(:);

%% stack dbn and softmax
stack = cell(numel(sizes),1);
for i=1:numel(sizes)
    stack{i}.w=dbn.rbm{i}.W;
    stack{i}.b=dbn.rbm{i}.c;
end
[stackparams, netconfig] = stack2params(stack);
stackedAETheta = [ saeSoftmaxOptTheta ; stackparams ];

%% Finetune the whole model
%opts.numepochs = 50;
stackedAEOptTheta=supervisedfinetuning(stackedAETheta, train_x, train_y,sizes(end),...
                numClasses, netconfig,softmaxlambda, opts,flag );


%% predict
acc=zeros(1,2);

[train_pred] = supervisedstackedAEPredict(stackedAETheta, sizes(end), numClasses, netconfig, train_x');
acc_temp = mean(train_y(:) == train_pred(:));
fprintf('Before Finetuning train Accuracy: %0.3f%%\n', acc_temp * 100);

[train_pred] = supervisedstackedAEPredict(stackedAEOptTheta, sizes(end), numClasses, netconfig, train_x');
acc_temp = mean(train_y(:) == train_pred(:));
fprintf('After Finetuning train Accuracy: %0.3f%%\n', acc_temp * 100);
acc(1)=acc_temp;

[test_pred] = supervisedstackedAEPredict(stackedAETheta, sizes(end), numClasses, netconfig, test_x');
acc_temp = mean(test_y(:) == test_pred(:));
fprintf('Before Finetuning test Accuracy: %0.3f%%\n', acc_temp * 100);

[test_pred] = supervisedstackedAEPredict(stackedAEOptTheta, sizes(end), numClasses, netconfig, test_x');
acc_temp = mean(test_y(:) == test_pred(:));
fprintf('After Finetuning test Accuracy: %0.3f%%\n', acc_temp * 100);
acc(2)=acc_temp;

output1=sizes(end);
output2=numClasses;
output3=netconfig;

% er=1-acc;

T=toc(ticid);
fprintf('The total time: %.4f\n',T);

end

