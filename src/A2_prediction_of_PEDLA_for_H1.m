function [ ] = A2_prediction_of_PEDLA_for_H1( neighbour, sizes, numepochs, cell_line )
%   Detailed explanation goes here   

%   This is the prediction procedure.
%   This function predicts the state of the unlabeled data based on the
%   learned model of PEDLA in the training procedure. All the input
%   paramerters should be same with the training procedure. 
%   The optimal paramerters we used in this paper
%   are  neighbour=0;sizes=[500 500];numepochs=50;cell_line='H1_cell_line'.
%   Since the training procedure will take about several hours and 40 G
%   physical memory, we use neighbour=0;sizes=[10];numepochs=5;cell_line='H1_cell_line' as default input parameters 
%   for quick investigation or test. You can omit the input parameter,
%   which will use the default input parameters.


%input:
%The unlabeled data is loaded from the default directory 
%'Data_for_training_and_prediction.H1/unlabelled_data-step200win200/H1_cell_line'
%and named 'chrN_neibour0.mat'.
%
%The trained model is loaded from the default directory
%'Model_learnt.H1/Model_learnt_q100/Model_learnt-ratio9step200win200/H1_cell_line'
%and named 'HMMSupervised_result_i50_500_500_n0.mat'. The learned model was trained in the training procedure.
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
%The prediction of the whole genome based on the trained model is saved in the default directory
%'Result_of_prediction.H1/H1_cell_line'
%and named 'PEDLA_prediction_of_H1_cell_line.chr1.txt'. Each of row
%is the state (label) of sequential 200-bp intevals on the corresponding
%chromosome with 1 indicating enhancer and 2 indicating non-enhancer.



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
my_quantile=100;  %normlization by this quantile.

disp(cell_line);


disp( ['numepochs=' num2str(opts.numepochs) '; sizes=[' num2str(sizes) ']; neighbour=' num2str(neighbour) '; ratio=' num2str(ratio) '; resolution=' num2str(resolution)  '; my_quantile=' num2str(my_quantile)]);

fprintf('\n');

dir_in=strcat('Model_learnt.H1');
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
dir_in=strcat(dir_in,'/',cell_line);
if ~exist(dir_in,'dir')
    mkdir(dir_in);
end

addpath(genpath('./Matlab_code_for_PEDLA/'));
dir_data='Data_for_training_and_prediction.H1';

%% loading quantile for normlization of unlabelled data
tic;
load(strcat(dir_data,'/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/q',num2str(my_quantile),'_of_',cell_line,'_neibour',num2str(neighbour),'.mat'),'q');   % load the quantile used to scale data
ttt1=toc;
disp(['loading the quantile takes time:' num2str(ttt1) 's']);


%% loading trained model of PEDLA
size_flag=[];
for i=1:length(sizes)
    size_flag=strcat(size_flag,'_',num2str(sizes(i)));
end
tic;
load(strcat(dir_in,'/HMMSupervised_result','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.mat'),'PI','A','B_DNN','T');    % load the learned model
ttt1=toc;
disp(['loading the trained model takes time:' num2str(ttt1) 's']);


%% output the predicted label of unlabelled data
prediction_result_dir='Result_of_prediction.H1/';
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