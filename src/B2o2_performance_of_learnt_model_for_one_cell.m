function [ performance_test ] = B2o2_performance_of_learnt_model_for_one_cell( neighbour,PI,A,B_DNN,cell_line )
%   Detailed explanation goes here

%   This function evaluates the performance of a trained PEDLA 
%   on one cell/tissue assigned by 'cell_line'.

%input:
%neighbour: the number of neighbours on both sides of the original one.
%       This is our inner testing parameter. 
%       In this paper, we set 'neighbour=0', and don't change it. 
%PI:    initial state probability vecort of HMM of a trained PEDLA 
%A: state transition probability matrix of HMM of a trained PEDLA 
%B_DNN: all model parameters of DNN of a trained PEDLA 
%cell_line: name of the cell line/tissue for evaluating.

%output:
%performance_test:    indicates the 6 performance indicators of a cell/tissue. 
%           Each performance indicator has two elements with
%           the former indicating the performance of the inner DNN and the
%           later indicating the performance of PEDLA. For
%           example, performance_test.acc(2) indicates the accuray of
%           PEDLA and performance_test.acc(1) indicates the accuray of the inner DNN. We don't
%           care the reuslt of the inner DNN here, so we need only view the
%           performance of PEDLA.


ratio=9;    %ratio of random regions/enhancers. alway set to be 9,except for the clss-imbanced analysis.
resolution=200; % resolution, which means the genome is divided into 200 bp intervals.
my_quantile=100;  %normlization by this quantile.

fprintf('\n');
disp(cell_line);

disp( ['ratio=' num2str(ratio) '; resolution=' num2str(resolution)  '; my_quantile=' num2str(my_quantile)]);

addpath(genpath('./Matlab_code_for_PEDLA/'));

%% loading data, pre-process and assignment of test set
tic;
load(strcat('Data_for_training_and_prediction/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/',cell_line,'_neibour',num2str(neighbour),'.mat'));
ttt1=toc;
disp(['loading takes time:' num2str(ttt1) 's']);

tic;
load(strcat('Data_for_training_and_prediction/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/q',num2str(my_quantile),'_of_twentytwo_training_cells_neibour',num2str(neighbour),'.mat'));
for kkk=1:numel(seqs)
    seqs{kkk} = bsxfun(@rdivide,seqs{kkk},q);    % normlization
end
ttt1=toc;
disp(['pre-processing takes time:' num2str(ttt1) 's']);

tic;
test_x=seqs;
test_y=states;
clear seqs;
clear states;
ttt1=toc;
disp(['assignment takes time:' num2str(ttt1) 's']);


%% performance of test set
num_algorithm=2;

TP=zeros(num_algorithm,1);
FP=zeros(num_algorithm,1);
TN=zeros(num_algorithm,1);
FN=zeros(num_algorithm,1);

tic;
real_states=[];
predicted_states=cell(num_algorithm,1);
for kk=1:numel(test_x)
    
    real_states=[real_states, test_y{kk}];
    
    [path,M2] = supervisedstackedAEPredict(B_DNN.stackedAEOptTheta,B_DNN.output1,B_DNN.output2,B_DNN.output3, test_x{kk});% inner DNN
    predicted_states{1}=[predicted_states{1},path];

    B = bsxfun(@rdivide,M2 ,B_DNN.p) ; %PEDLA
    [path] = lf_HMMViterbi(PI,A, B);
    predicted_states{2}=[predicted_states{2},path];
end
ttt2=toc;
disp(['prediction of test set takes time:' num2str(ttt2) 's']);

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
performance_test.acc=acc;
performance_test.sen=sen;
performance_test.spe=spe;
performance_test.gm=gm;
performance_test.ppv=ppv;
performance_test.f1=f1;
disp(strcat('test accuracy of DNN:',num2str(acc(1)*100),'%')) 
disp(strcat('test accuracy of PEDLA:',num2str(acc(2)*100),'%'))





end

