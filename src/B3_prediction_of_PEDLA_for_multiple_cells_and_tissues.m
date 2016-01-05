function [ ] = B3_prediction_of_PEDLA_for_multiple_cells_and_tissues( neighbour, sizes, numepochs, order_ind, N )
%   Detailed explanation goes here   

%   This is the prediction procedure. This function is almost the same with
%   'A2_prediction_of_PEDLA_for_H1.m' and uses the model trained on multiple cells/tissues.
%   This function predicts the state of the unlabeled data based on the
%   learned model of PEDLA trained on 22 cells/tissues in the training procedure.
%
%   The optimal paramerters we used for training the PEDLA in our paper are neighbour=0;sizes=[50];numepochs=150;order_ind=1;N=22.
%   And the input paramerters should be consistent with those of the
%   training procedure. Since the training procedure is extremly time-consuming
%   and finishing training on the whole 22 training cells/tissues may take a few days, 
%   we use neighbour=0;sizes=[5];numepochs=5;order_ind=1;N=22 as default input parameters 
%   for quick investigation or test, which can finish training qucikly. You can omit the input parameter,
%   which will use the default input parameters.


%input:
%The unlabeled data is loaded from sub-directories of the default directory 
%'Data_for_training_and_prediction/unlabelled_data-step200win200/'
%and named 'chrM_neibour0.mat'. M is chromosome index.
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


%output
%The prediction of the whole genome of the corresponding cells/tissues based on the trained model
%is saved in sub-directories of the default directory 'Result_of_prediction'
%and named 'PEDLA_prediction_of_XXX.chr1.txt'. XXX indicates the name of cell/tissue.
%Each row of the result is the state (label) of sequential 200-bp intevals on the corresponding
%chromosome with 1 indicating enhancer and 2 indicating non-enhancer.



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


ratio=9;    %ratio of random regions/enhancers. alway set to be 9,except for the clss-imbanced analysis.
resolution=200; % resolution, which means the genome is divided into 200 bp intervals.
my_quantile=100;  %normlization by this quantile.

disp(['order' num2str(order_ind) ' and model ' num2str(N)]);


disp( ['numepochs=' num2str(opts.numepochs) '; sizes=[' num2str(sizes) ']; neighbour=' num2str(neighbour) '; ratio=' num2str(ratio) '; resolution=' num2str(resolution)  '; my_quantile=' num2str(my_quantile)]);

fprintf('\n');

dir_in=strcat('Model_learnt');
if ~exist(dir_in,'dir')
    mkdir(dir_in);
end
dir_in=strcat(dir_in,'/Model_learnt_q',num2str(my_quantile));
if ~exist(dir_in,'dir')
    mkdir(dir_in);
end
dir_in=strcat(dir_in,'/Model_learnt-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution));
if ~exist(dir_in,'dir')
    mkdir(dir_in);
end
dir_in=strcat(dir_in,'/Order',num2str(order_ind));
if ~exist(dir_in,'dir')
    mkdir(dir_in);
end

addpath(genpath('./Matlab_code_for_PEDLA/'));
dir_data='Data_for_training_and_prediction';

%% loading quantile for normlization of unlabelled data
tic;
load(strcat(dir_data,'/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/q',num2str(my_quantile),'_of_twentytwo_training_cells_neibour',num2str(neighbour),'.mat'),'q');   % load the quantile used to scale data
ttt1=toc;
disp(['loading the quantile takes time:' num2str(ttt1) 's']);


%% loading trained model of PEDLA
size_flag=[];
for i=1:length(sizes)
    size_flag=strcat(size_flag,'_',num2str(sizes(i)));
end
tic;
load(strcat(dir_in,'/HMMSupervised_result','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.',num2str(N),'.mat'),'PI','A','B_DNN','T'); % load the trained PEDLA 
ttt1=toc;
disp(['loading the trained model takes time:' num2str(ttt1) 's']);


%% output the predicted label of unlabelled data
prediction_result_dir='Result_of_prediction/';
if ~exist(prediction_result_dir,'dir')
    mkdir(prediction_result_dir);
end
unlabelled_dir=strcat(dir_data,'/unlabelled_data-','step',num2str(resolution),'win',num2str(resolution),'/');
unlabelled_cells = dir(unlabelled_dir);
for ii=3:numel(unlabelled_cells)
    disp(['writing results of ' unlabelled_cells(ii).name '...']);
    unlabelled_cell = unlabelled_cells(ii).name;
    
    if ~exist(strcat(prediction_result_dir,unlabelled_cell),'dir')
        mkdir(strcat(prediction_result_dir,unlabelled_cell));
    end
    
    %%load unlabelled data for prediction
    for i=1:23
        chr=strcat('chr',num2str(i));
        if i==23
            chr='chrX';
        end
        disp( ['  ',chr] );
        
        
        load(strcat(unlabelled_dir,unlabelled_cell,'/',chr,'_neibour0.mat'),'unlabelled_seqs'); % load the unlabeled data
        unlabelled_seqs = bsxfun(@rdivide,unlabelled_seqs,q);   % normlization
        
        [~,M2] = supervisedstackedAEPredict(B_DNN.stackedAEOptTheta,B_DNN.output1,B_DNN.output2,B_DNN.output3, unlabelled_seqs);%output of DNN
        B = bsxfun(@rdivide,M2 ,B_DNN.p) ; 
        [path] = lf_HMMViterbi(PI,A, B);    %prediction of PEDLA
        
        disp('  writing results of PEDLA...')
        dlmwrite(strcat(prediction_result_dir,unlabelled_cell,'/PEDLA_prediction_of_',unlabelled_cell,'.',chr,'.txt'),path(:),'delimiter','' );

    end
end



end