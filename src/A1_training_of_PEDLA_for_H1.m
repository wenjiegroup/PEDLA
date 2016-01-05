function [ T, record ] = A1_training_of_PEDLA_for_H1( neighbour, sizes, numepochs, cell_line )
%   Detailed explanation goes here

%   This is the training procedure.
%   This function learns the model parameters of PEDLA. The optimal paramerters we used in this paper
%   are  neighbour=0;sizes=[500 500];numepochs=50;cell_line='H1_cell_line'.
%   Since the training procedure will take about several hours and 40 G
%   physical memory, we use neighbour=0;sizes=[10];numepochs=5;cell_line='H1_cell_line' as default input parameters 
%   for quick investigation or test. You can omit the input parameter,
%   which will use the default input parameters.

%   The training procedure will not get a same result with the same input
%   parameters because it used a different random seed related with the current
%   time each time we run it. To get a same result with the same input
%   parameters, you can change the code below from rng('shuffle') to rng(0).


%input:
%The input data is loaded from the default directory 
%'Data_for_training_and_prediction.H1/Data_for_training_and_prediction-ratio9step200win200'
%and named 'H1_cell_line_neibour0.mat'. You can change it by assigh the
%input parameter cell_line, since cell_line is default 'H1_cell_line'.  
%
%neighbour: the number of neighbours on both sides of the original one.
%       This is our inner testing parameter. 
%       In this paper, we set 'neighbour=0', and don't change it. 
%sizes: structure of PEDLA. A vector indicating the numbers of units in each
%       hidden layer. In this paper ,we use a optimal structure of [500 500] for H1 cell line with 1114-dimensional features,
%       which means PEDLA has 2 hidden layers and each layer has 500 units.
%numepochs: number of epochs for PEDLA. In this paper ,we use 50 here.
%cell_line: name of the cell line/tissue for training. Here is
%       'H1_cell_line' for our paper. You can change it to other
%       cell/tissue, if you have corresponding input data.


%output

%The trained model is saved in the default directory
%'Model_learnt.H1/Model_learnt_q100/Model_learnt-ratio9step200win200/H1_cell_line'
%and named 'HMMSupervised_result_i50_500_500_n0.mat'. The learned model
%can be used for later prediction and analysis.
%
%T:     time cost in second. This is additional information.
%record:    record of performance indicators of the training data for inner DNN and PEDLA. This is additional information.


if nargin<1
    neighbour = 0;   %the number of neibours on both sides of the original one. constant: 0
end
if nargin<2
    sizes = [10];   %structure of PEDLA. defualt: [10]; optimal: [500 500]
end
if nargin<3
    numepochs = 5;   % defualt: 5; optimal: 50
end
opts.numepochs=numepochs;
if nargin<4
	cell_line='H1_cell_line';    % defualt: 'H1_cell_line'. You can change it to other cell/tissue, if you have corresponding input data
end


ratio=9;    %ratio of random regions/enhancers. alway set to be 9,except for the clss-imbanced analysis.
resolution=200; % resolution, which means the genome is divided into 200 bp intervals.
my_quantile=100;  %scaled to [0,1] by this quantile.

disp(cell_line);


disp( ['numepochs=' num2str(opts.numepochs) '; sizes=[' num2str(sizes) ']; neighbour=' num2str(neighbour) '; ratio=' num2str(ratio) '; resolution=' num2str(resolution)  '; my_quantile=' num2str(my_quantile)]);
fprintf('\n');

dir_out=strcat('Model_learnt.H1');
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end
dir_out=strcat(dir_out,'/Model_learnt_q',num2str(my_quantile));
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end
dir_out=strcat(dir_out,'/Model_learnt-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution));
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end
dir_out=strcat(dir_out,'/',cell_line);
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end

addpath(genpath('./Matlab_code_for_PEDLA/'));


%% loading input data, pre-processing and assignment of trainging set
tic;
load(strcat('Data_for_training_and_prediction.H1/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/',cell_line,'_neibour',num2str(neighbour),'.mat'),'seqs','states');
ttt1=toc;
disp(['loading takes time:' num2str(ttt1) 's']);

tic;
temp=[];
for k = 1:numel(seqs) 
    temp=[temp,seqs{k}];
end
q=quantile(temp',[my_quantile/100])';
save(strcat('Data_for_training_and_prediction.H1/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/q',num2str(my_quantile),'_of_',cell_line,'_neibour',num2str(neighbour),'.mat'),'q');   % save the quantile used to scale data to [0,1]
temp=[];
for kkk=1:numel(seqs)
    seqs{kkk} = bsxfun(@rdivide,seqs{kkk},q);    % scaled to [0,1]
end
ttt1=toc;
disp(['pre-processing takes time:' num2str(ttt1) 's']);

tic;
training_x=seqs;
training_y=states;
clear seqs;
clear states;
ttt1=toc;
disp(['assignment takes time:' num2str(ttt1) 's']);


%% learning the whole model of PEDLA here
rng('shuffle'); % we set 'shuffle' to produce a different sequence of random numbers based on the current time. set 0 to get a same random seed.

[ PI,A,B_DNN,T,acc_DNN ] = lf_HMMSupervisedLearning( training_x, training_y, training_x, training_y, opts, sizes);  % the main learning function

size_flag=[];
for i=1:length(sizes)
    size_flag=strcat(size_flag,'_',num2str(sizes(i)));
end
save(strcat(dir_out,'/HMMSupervised_result','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.mat'),'PI','A','B_DNN','T');   % save the learned model


%% performance evaluation of the learned model on the training set
% The codes below are optional, since the above codes have learned the
% model parameters of PEDLA. We assess the performance of the learned model
% on the same training data to get a glimpse of the performance.
num_algorithm=2;

TP=zeros(num_algorithm,1);
FP=zeros(num_algorithm,1);
TN=zeros(num_algorithm,1);
FN=zeros(num_algorithm,1);

tic;
real_states=[];
predicted_states=cell(num_algorithm,1);
for kk=1:numel(training_x)
    
    real_states=[real_states, training_y{kk}];
  
    [path,M2] = supervisedstackedAEPredict(B_DNN.stackedAEOptTheta,B_DNN.output1,B_DNN.output2,B_DNN.output3, training_x{kk});% inner DNN
    predicted_states{1}=[predicted_states{1},path];

    B = bsxfun(@rdivide,M2 ,B_DNN.p) ; %PEDLA
    [path] = lf_HMMViterbi(PI,A, B);
    predicted_states{2}=[predicted_states{2},path];
end
ttt2=toc;
disp(['prediction of training set takes time:' num2str(ttt2) 's']);

if size(predicted_states{1},2)~=size(predicted_states{2},2) || size(predicted_states{1},2)~=size(real_states,2)
    error('size not equal!');
end
for i=1:numel(predicted_states)
    TP(i)=TP(i)+length(find(real_states==1 & predicted_states{i}==1));
    FP(i)=FP(i)+length(find(real_states==2 & predicted_states{i}==1));
    TN(i)=TN(i)+length(find(real_states==2 & predicted_states{i}==2));
    FN(i)=FN(i)+length(find(real_states==1 & predicted_states{i}==2));
end
acc=(TP+TN)./(TP+TN+FP+FN);
sen=TP./(TP+FN);
spe=TN./(TN+FP);
gm=sqrt(sen.*spe);
ppv=TP./(TP+FP);
f1=2./(1./sen+1./ppv);
performance_training.acc=acc;
performance_training.sen=sen;
performance_training.spe=spe;
performance_training.gm=gm;
performance_training.ppv=ppv;
performance_training.f1=f1;
disp(strcat('training accuracy of DNN:',num2str(acc(1)*100),'%'))
disp(strcat('training accuracy of PEDLA:',num2str(acc(2)*100),'%'))


record=[performance_training];

save(strcat(dir_out,'/performance','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.mat'),'record');     % save performance of the learned model on the same training set





end