function [  ] = B1_training_of_PEDLA_for_multiple_cells_and_tissues( neighbour, sizes, numepochs, order_ind )
%   Detailed explanation goes here

%   This is the training procedure of PEDLA for multiple cells and tissues.
%   This function learns the model parameters of PEDLA from multiple cells and tissues. The
%   PEDLA for multiple cells and tissues is a general framework of PEDLA and can be used to train and predicte
%   on multiple cell lines/tissues. PEDLA is a special case of
%   PEDLA for multiple cells and tissues where the number of training cells/tissues used for PEDLA is 1.
%   The optimal paramerters we used in this paper are  neighbour=0;sizes=[50];numepochs=150;order_ind=1.
%   Since the training procedure is extremly time-consuming and finishing training on the whole 22 training cells
%   and tissues may take a few days, we use neighbour=0;sizes=[5];numepochs=5;order_ind=1 as default input parameters 
%   for quick investigation or test, which can finish training very qucikly. You can omit the input parameter,
%   which will use the default input parameters.

%   The training procedure will not get a same result with the same input
%   parameters because it used a different random seed related with the current
%   time each time we run it. To get a same result with the same input
%   parameters, you can change the code below from rng('shuffle') to rng(0).


%input:
%The input data is loaded from the default directory 
%'Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200'
%by the order of training cells/tissues. And the order of training
%cells/tissues is loaded from the default directory
%'Order_of_22_training_cells_and_tissues' and named 'Order1.txt'.
%
%neighbour: the number of neighbours on both sides of the original one.
%       This is our inner testing parameter. 
%       In this paper, we set 'neighbour=0', and don't change it. 
%sizes: structure of PEDLA. A vector indicating the numbers of units in each
%       hidden layer. In our paper ,we use a optimal structure of [50] for training 
%       PEDLA on 22 training cells and tissues,
%       which means PEDLA has 1 hidden layer and the layer has 50 units.
%numepochs: number of epochs for PEDLA. In our paper ,we use 150 here.
%order_ind: the index of order of the training cells/tissues. defualt and optimal: 1


%output

%The trained model is saved in the default directory
%'Model_learnt/Model_learnt_q100/Model_learnt-ratio9step200win200/Order1'
%and named 'HMMSupervised_result_i150_50_n0.N.mat'. N is the number of
%cells/tissues that PEDLA have trained on. 
%The learned model can be used for later training, prediction and evaluation.
%



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



ratio=9;    %ratio of random regions/enhancers. alway set to be 9,except for the clss-imbanced analysis.
resolution=200; % resolution, which means the genome is divided into 200 bp intervals.
my_quantile=100;  %normlization by this quantile.

fileID = fopen( strcat('Order_of_22_training_cells_and_tissues/Order',num2str(order_ind),'.txt'), 'r' );
order_cell=textscan(fileID,'%d\t%s');   % read the order of training cells/tissues.
fclose(fileID);
number_of_cell=numel(order_cell{2}); 

disp(['order' num2str(order_ind)]);

for cell_ind=1:number_of_cell
    
    cell_line=order_cell{2}{cell_ind};
    fprintf('\n');
    disp(cell_line);

    disp( ['numepochs=' num2str(opts.numepochs) '; sizes=[' num2str(sizes) ']; neighbour=' num2str(neighbour) '; ratio=' num2str(ratio) '; resolution=' num2str(resolution)  '; my_quantile=' num2str(my_quantile)]);

    dir_out=strcat('Model_learnt');
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
    dir_out=strcat(dir_out,'/Order',num2str(order_ind));
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end

    addpath(genpath('./Matlab_code_for_PEDLA/'));

    %% loading data, pre-processing and assignment of trainging set
    tic;
    load(strcat('Data_for_training_and_prediction/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/',cell_line,'_neibour',num2str(neighbour),'.mat'),'seqs','states');
    ttt1=toc;
    disp(['loading takes time:' num2str(ttt1) 's']);

    tic;
    load(strcat('Data_for_training_and_prediction/Data_for_training_and_prediction-ratio',num2str(ratio),'step',num2str(resolution),'win',num2str(resolution),'/q',num2str(my_quantile),'_of_twentytwo_training_cells_neibour',num2str(neighbour),'.mat'),'q');
    for kkk=1:numel(seqs)
        seqs{kkk} = bsxfun(@rdivide,seqs{kkk},q);    % normlization
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

    %% learning the whole model of PEDLA for multiple cells and tissues here
    rng('shuffle'); % we set 'shuffle' to produce a different sequence of random numbers based on the current time. set 0 to get a same random seed.
    size_flag=[];
    for i=1:length(sizes)
        size_flag=strcat(size_flag,'_',num2str(sizes(i)));
    end

    if cell_ind==1
        [ PI,A,B_DNN,T,acc_DNN ] = lf_HMMSupervisedLearning( training_x, training_y, training_x, training_y, opts, sizes);  % the initial training of PEDLA for multiple cells and tissues
    elseif cell_ind>1
        load(strcat(dir_out,'/HMMSupervised_result','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.',num2str(cell_ind-1),'.mat'),'PI','A','B_DNN','T');   % fetch the trained model parameters of a previous cell line/tissue as initialization
        [ PI,A,B_DNN,T,acc_DNN ] = lf_HMMSupervisedLearning_FromExistedModel( training_x, training_y, training_x, training_y, opts, sizes, PI, A, B_DNN );  % the iterative training of PEDLA for multiple cells and tissues
    end

    save(strcat(dir_out,'/HMMSupervised_result','_i',num2str(opts.numepochs),size_flag,'_n',num2str(neighbour),'.',num2str(cell_ind),'.mat'),'PI','A','B_DNN','T'); % save the learned model

end



end