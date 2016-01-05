function [ PI_new, A_new, B_DNN_new, T, acc_DNN ] = lf_HMMSupervisedLearning_FromExistedModel( training_x, training_y, test_x, test_y, opts, sizes, PI_old, A_old, B_DNN_old )
%LF_LEARNINGFROMEXISTEDMODEL Summary of this function goes here
%   Detailed explanation goes here


if ~iscell(training_x) || ~iscell(training_y)
    error(message('stats:lf_HMMSupervisedLearning:Type of both inputs must be cell'));    
end

K = numel(training_x);
if K ~= numel(training_y)
    error(message('stats:lf_HMMSupervisedLearning:Number of input sequences must equal'));
end

for i=1:K
    if size(training_x{i},2) ~= size(training_y{i},2)
        error(message('stats:lf_HMMSupervisedLearning:Length of observation and state must equal'));
    end
end
   
D = size(training_x{1},1);
for i=1:K
    if size(training_x{i},1) ~= D
        error(message('stats:lf_HMMSupervisedLearning:Dimension of observation must equal'));
    end
end

if nargin < 5
    opts.numepochs = 1;
end
if nargin < 6
    sizes = [200 200];
end



%% initiate pararmeter
disp('initiating pararmeter...');
AllState = [];
for i=1:K    
    AllState = [AllState, unique(training_y{i})];
end
N = length(unique(AllState));
PI = zeros(N,1);
A = zeros(N);
% B1.info = 'GMM';
% B1.mixmat = zeros(N,M);
% B1.mu = zeros(D,N,M);
% B1.sigma = zeros(D,D,N,M);
% B3.info = 'GMM';

pseudoPI  = ones(size(PI))-1;
pseudoA = ones(size(A))-1;

%% calculate PI and A
disp('calculating PI and preparing data for next step...');
training_data=[];
training_label=[];

tic;
for i = 1:N
    data{i} = [];
end
for k = 1:K 
%     if mod(k,1000)==1
%         disp(strcat('  sequence:',num2str(k)));
%     end
%     Tk = length(training_y{k});
    PI(training_y{k}(1)) = PI(training_y{k}(1)) + 1;     
%     for t=1:Tk-1
%         A(training_y{k}(t),training_y{k}(t+1)) = A(training_y{k}(t),training_y{k}(t+1)) + 1;   
%         % xi{k}(:,:,t) =      % full( sparse(training_y{k}(1:Tk-1) ,training_y{k}(2:Tk) ,1,N,N) );
%         % gamma{k}(:,t) = 
%     end
%     for i = unique(training_y{k})
%         data{i} = [data{i},training_x{k}(:,find(training_y{k}==i))];
%     end
    training_data=[training_data,training_x{k}];
    training_label=[training_label,training_y{k}];
end
for i = unique(training_label)
    data{i} = [data{i},training_data(:,find(training_label==i))];
end;
ttt1=toc;

clear training_x training_y;
disp(['calculating PI and preparing data takes time:' num2str(ttt1) 's']);
fprintf('\n');

PI = PI + pseudoPI;
A = A + pseudoA;

PI = PI/sum(PI);
% A = bsxfun(@rdivide,A,sum(A,2));
A=[];%训练集的特殊性导致学习到的A必然是[1 0; 0 1],不可能学习到真正的A，所以这里我们就不计算A。因为计算A的外层循环是10^4级别，很耗时，至少几个小时都没跑完；即使计算了也没用，因为计算的是错的。

PI_new = PI_old + PI;
PI_new = PI_new/sum(PI_new);
A_new=[];

T=[];

%% DNN
% addpath('.');
% [~, model] = emgm(training_data, M);   % there is no need to calculate P(Xt), just for checking. In final version, this can be omitted
% obj = gmdistribution(model.mu',model.Sigma,model.weight); 

test_data=[];
test_label=[];
for k = 1:numel(test_y)   %%  get test data and label
    test_data=[test_data,test_x{k}];
    test_label=[test_label,test_y{k}];
end

tic;
[ acc_DNN,stackedAEOptTheta,output1,output2,output3,T_DNN ] = supervisedlearnig_from_existed_model(training_data',training_label',test_data',test_label',sizes, 0, opts,B_DNN_old); %% the main learning function of DNN
T(end+1)=T_DNN;
B_DNN_new.stackedAEOptTheta = stackedAEOptTheta;
B_DNN_new.output1 = output1;
B_DNN_new.output2 = output2;
B_DNN_new.output3 = output3;
for i=unique(AllState)
    B_DNN.p(i) = (size(data{i},2)/size(training_data,2));
end
B_DNN.p =  B_DNN.p(:);
B_DNN_new.p = B_DNN_old.p + B_DNN.p;
B_DNN_new.p = B_DNN_new.p/sum(B_DNN_new.p);
% B_DNN.obj = obj;
B_DNN_new.info = 'DNN';
T(end+1)=toc;


end



function [ acc,stackedAEOptTheta,output1,output2,output3, T ] = supervisedlearnig_from_existed_model(train_x,train_y,test_x,test_y, sizes, flag, opts,B_DNN_old )

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

% %% unsupervised feature learning
% dbn = dbnsetup(sizes, train_x, opts);
% dbn = dbntrain(dbn, train_x, opts);
% 
% %% input of softmax
% trainfeatue=train_x;
% for i=1:numel(sizes)
%     trainfeatue=rbmup(dbn.rbm{i},trainfeatue);
% end
% inputSize=size(trainfeatue,2);
% numClasses=length(unique(train_y));
% softmaxModel=softmaxtrain(inputSize,numClasses , trainfeatue, train_y, opts, softmaxlambda,flag );
% %softmaxoutput=softmaxtup( softmaxModel.optTheta , trainfeatue');
% saeSoftmaxOptTheta = softmaxModel.optTheta(:);
% 
% %% stack dbn and softmax
% stack = cell(numel(sizes),1);
% for i=1:numel(sizes)
%     stack{i}.w=dbn.rbm{i}.W;
%     stack{i}.b=dbn.rbm{i}.c;
% end
% [stackparams, netconfig] = stack2params(stack);
% stackedAETheta = [ saeSoftmaxOptTheta ; stackparams ];

%% Finetune the whole model
%opts.numepochs = 50;

numClasses=B_DNN_old.output2;
netconfig=B_DNN_old.output3;
stackedAEOptTheta=supervisedfinetuning(B_DNN_old.stackedAEOptTheta, train_x, train_y,sizes(end),...
                numClasses, netconfig,softmaxlambda, opts,flag );


%% predict
acc=zeros(1,2);

[train_pred] = supervisedstackedAEPredict(B_DNN_old.stackedAEOptTheta, sizes(end), numClasses, netconfig, train_x');
acc_temp = mean(train_y(:) == train_pred(:));
fprintf('Before Finetuning train Accuracy: %0.3f%%\n', acc_temp * 100);

[train_pred] = supervisedstackedAEPredict(stackedAEOptTheta, sizes(end), numClasses, netconfig, train_x');
acc_temp = mean(train_y(:) == train_pred(:));
fprintf('After Finetuning train Accuracy: %0.3f%%\n', acc_temp * 100);
acc(1)=acc_temp;

[test_pred] = supervisedstackedAEPredict(B_DNN_old.stackedAEOptTheta, sizes(end), numClasses, netconfig, test_x');
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


