function [ PI_new, A_new, B_DNN_new, T, acc_DNN ] = lf_HMMSupervisedLearning_FromExistedModel( training_x, training_y, test_x, test_y, opts, sizes, PI_old, A_old, B_DNN_old )
%LF_LEARNINGFROMEXISTEDMODEL Summary of this function goes here
%   Detailed explanation goes here
%   the main learning process of PEDLA. DNN is a part of PEDLA.

%
%input:
%training_x:  the training feature. A cell vector with each cell indicating
%               the features of an ehnhancer or a non-enhancer. For example, if
%               "chr1 1000 2200" is enhancer, then its corresponding cell
%               is a D*6 matrix. D is the dimension of the feature(for
%               exmaple 1114).Six is the lenght of "chr1 1000 2200", since
%               (2200-1000)/200=6. 200 is the resolution,which means the genome is divided into 200 bp intervals.
%training_y:  the training lable. Similar with training_x. A cell vector with each cell indicating
%               the labels of an ehnhancer or a non-enhancer (1 for enhancer, 2 for non-enhancer). For example, if
%               "chr1 1000 2200" is enhancer, then its corresponding cell
%               is a [1 1 1 1 1 1] vector.
%test_x:    the test feature. same structure as training_x. training_x and
%           training_y are used to train PEDLA; test_x and test_y are used to
%           evaluate PEDLA.
%test_y:    the test lable. same structure as training_y.
%opts:  other optional parameter.
%size:  structure of PEDLA
%PI_old:    the initial state probability vecort of a previous model.
%A_old:         state transition probability matrix of a previous model.
%B_DNN_old:     all model parameters of a previous model.
%


%output:
%PI_new:    refined initial state probability vecort of HMM
%A_new: refined state transition probability matrix of HMM.
%B_DNN_new: all refined model parameters of DNN
%T: time cost of these algorithms
%acc_DNN:   traning and test accuracy of DNN

% addpath(genpath('./'));

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
N = length(unique(AllState));   % number of all states. For enhancer prediction, N is always 2.
PI = zeros(N,1);
A = zeros(N);


pseudoPI  = ones(size(PI))-1;
pseudoA = ones(size(A))-1;

%% calculate PI and A of HMM
disp('calculating PI and A, and preparing data for next step...');
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
%     end
    training_data=[training_data,training_x{k}];
    training_label=[training_label,training_y{k}];
end
for i = unique(training_label)
    data{i} = [data{i},training_data(:,find(training_label==i))];
end;
ttt1=toc;

clear training_x training_y;
disp(['calculating PI A and preparing data takes time:' num2str(ttt1) 's']);
fprintf('\n');

PI = PI + pseudoPI;
A = A + pseudoA;

PI = PI/sum(PI);
% A = bsxfun(@rdivide,A,sum(A,2));

%NOTE: As shown in Fig. S1, the algorithm automatically learns the triple (PI, A, B) in the training procedure. 
%However, due to the discreteness of training data the state transition probability matrix A cannot be estimated correctly (the code is right, and if the training data surports and the A will be estimated correctly), 
%where A would constantly be matrix [1 0; 0 1] if learning automatically. 
%The discreteness of training data refers to that there is not a continuous genomic locus in the training data where two states (enhancer class and non-enhancer class) transit explicitly, 
%which is caused by the fact that there is no clear boundary between
%enhancer and non-enhancer states. If the training data are like this,
%"chr1 1000 2200" is enhancer, "chr1 2200 5000" is non-enhancer, "chr1 5000
%7000" is enhancer..., the code can lean the A correctly.
%We commented out the codes above that is used to learn A, since it cost time and the learned A can not reflect
%the real situation.
%For the mentioned reason, we had to designate A manually.A was constantly set to [0.95 0.05; 0.0005 0.9995] empirically throughout this paper, 
%which means the probability of transition from enhancer class to enhancer class, transition from enhancer class to non-enhancer class, 
%transition from non-enhancer class to enhancer class and transition from non-enhancer class to non-enhancer class is 0.95, 0,05, 0.0005 and 0.9995 respectively.
A=[19  1 ; 1  1999  ];
A = bsxfun(@rdivide,A,sum(A,2));

PI_new = PI_old + PI;
PI_new = PI_new/sum(PI_new);
A_new=  A_old+A;
A_new = bsxfun(@rdivide,A_new,sum(A_new,2));

T=[];

%% DNN
test_data=[];
test_label=[];
for k = 1:numel(test_y)   %%  get test data and label
    test_data=[test_data,test_x{k}];
    test_label=[test_label,test_y{k}];
end
clear test_x test_y;

tic;
[ acc_DNN,stackedAEOptTheta,output1,output2,output3,T_DNN ] = supervisedlearnig_from_existed_model(training_data',training_label',test_data',test_label',sizes, 0, opts,B_DNN_old);%% the main learning function of DNN
T(end+1)=T_DNN;
B_DNN_new.stackedAEOptTheta = stackedAEOptTheta;     % the weight matirx W of DNN
B_DNN_new.output1 = output1;
B_DNN_new.output2 = output2;
B_DNN_new.output3 = output3;
for i=unique(AllState)
    B_DNN.p(i) = (size(data{i},2)/size(training_data,2));    % the prior probability P of PEDLA used for alleviating class-imbalance issue
end
B_DNN.p =  B_DNN.p(:);
B_DNN_new.p = B_DNN_old.p + B_DNN.p;
B_DNN_new.p = B_DNN_new.p/sum(B_DNN_new.p);
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


