function [  ] = B2_performance_evaluation_of_on_training_or_test_set( neighbour, sizes, numepochs, order_ind, N, test_or_training )
%   Detailed explanation goes here

%   This is performance evaluation of the PEDLA on our 22 training cells/tissues or 20 independent test cells/tissues.
%   This function evaluates the performance of a trained PEDLA and is optional, since the training and
%   prediction procedure of the PEDLA are necessary.
%   
%   The optimal paramerters we used for training the PEDLA for multiple cells and tissues in our paper are neighbour=0;sizes=[50];numepochs=150;order_ind=1;N=22.
%   And the input paramerters should be consistent with those of the
%   training procedure. Since the performance evaluation is time-consuming, 
%   we use neighbour=0;sizes=[5];numepochs=5;order_ind=1;N=22 as default input parameters 
%   for quick investigation or test, which can finish evaluating qucikly. You can omit the input parameter,
%   which will use the default input parameters.


%input:
%The input data for evaluation is loaded from the default directory 
%'Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200'
%
%The trained model is loaded from the default directory
%'Model_learnt/Model_learnt_q100/Model_learnt-ratio9step200win200/Order1'
%and named 'HMMSupervised_result_i150_50_n0.N.mat'. N is the number of
%cells/tissues that PEDLA have trained on. 
%
%neighbour: the number of neighbours on both sides of the original one.
%       This is our inner testing parameter. 
%       In this paper, we set 'neighbour=0', and don't change it. 
%sizes: structure of PEDLA. A vector indicating the numbers of units in each
%       hidden layer. In our paper ,we use a optimal structure of [50] for training 
%       PEDLA on 22 training cells and tissues,
%       which means PEDLA has 1 hidden layer and the layer has 50 units.
%numepochs: number of epochs for PEDLA. In this paper ,we use 150 here.
%order_ind: the index of order of the training cells/tissues. defualt and optimal: 1. 
%           should be consistent with that of the training procedure
%N: number of cells/tissues that PEDLA have been trained on. can be 1~22. defualt and optimal: 22
%test_or_training:  should be either 'test' or 'training' which indicates evaluation on test or training cells/tissues.
%                   defualt: 'test'


%output

%The result of performance evaluation is saved in the default directory
%'Model_learnt/Model_learnt_q100/Model_learnt-ratio9step200win200/Order1'
%and named 'performance_of_model_i150_50_n0.training_cells.N.mat' or 
%named 'performance_of_model_i150_50_n0.test_cells.N.mat'. 'N' is the number of
%cells/tissues that PEDLA have trained on. 'training_cells'
%indicates that the performance is evaluated on the 22 training cells/tissues 
%and 'test_cells' indicates that the performance is evaluated on 20 independent test cells/tissues.
%The result of performance evaluation is saved in a variable 'record' in
%file 'performance_of_model_i150_50_n0.training_cells.N.mat' or
%'performance_of_model_i150_50_n0.test_cells.N.mat'.
%record:    record of performance evaluation. A vector with each element
%           indicating the 6 performance indicators of a cell/tissue. The
%           length of the record equals the number of cells/tissues
%           evaluated on. Each performance indicator has two elements with
%           the former indicating the performance of the inner DNN and the
%           later indicating the performance of PEDLA. For
%           example, record(15).acc(2) indicates the accuray of
%           PEDLA on the 15-th cell/tissue and record(15).acc(1) indicates
%           the accuray of the inner DNN on the 15-th cell/tissue. We don't
%           care the reuslt of the inner DNN here, so we need only view the
%           performance of PEDLA.
%
%To evaluate the performance of the 22 sequentially trained model on 20 independent test cells/tissues, you can
%run the code "for i=1:22 B2_performance_evaluation_of_on_training_or_test_set( 0, [50], 150, 1, i, 'test' ); end"


if nargin<1
    neighbour = 0;   %the number of neibours on both sides of the original one. constant: 0
end
if nargin<2
    sizes = [5];   %structure of PEDLA. defualt: [5]; optimal: [50]
end
if nargin<3
    numepochs = 5;   % defualt: 5; optimal: 150
end
opts.numepochs=numepochs;
if nargin<4
	order_ind=1;    % the index of order of the training cells/tissues. defualt and optimal: 1
end
if nargin<5
	N=22;    % number of cells/tissues that PEDLA have trained on. can be 1~22. defualt and optimal: 22
end
if nargin<6
	test_or_training='test';    % test or training cells/tissues.
end


ratio=9;    %ratio of random regions/enhancers. alway set to be 9,except for the clss-imbanced analysis.
resolution=200; % resolution, which means the genome is divided into 200 bp intervals.
my_quantile=100;  %normlization by this quantile.

disp(['order' num2str(order_ind) ' and model ' num2str(N)]);

dir_in=strcat('Model_learnt');
dir_in=strcat(dir_in,'/Model_learnt_q',num2str(my_quantile));
dir_in=strcat(dir_in,'/Model_learnt-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution));
dir_in=strcat(dir_in,'/Order',num2str(order_ind));


size_flag=[];
for i=1:length(sizes)
    size_flag=strcat(size_flag,'_',num2str(sizes(i)));
end

load(strcat(dir_in,'/HMMSupervised_result','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.',num2str(N),'.mat'),'PI','A','B_DNN','T'); % load the trained PEDLA 

[ record ] = B2o1_merge_performance_of_training_or_test_set( PI,A,B_DNN,test_or_training );
save( strcat(dir_in,'/performance_of_model_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.',test_or_training,'_cells.',num2str(N),'.mat') , 'record' ); % save the performance



end

