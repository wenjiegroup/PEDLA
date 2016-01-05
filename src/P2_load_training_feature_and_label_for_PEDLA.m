function [ seqs,states ] = P2_load_training_feature_and_label_for_PEDLA( neibour, cell )
%   Detailed explanation goes here

%   This function just loads the extracted features and labels into a
%   matlab format.


%input:
%The features and labels extracted in a previous step are loaded from sub-directories of the default directory 
%'Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200/'.
%The files prefixed by "x_" are the features and the other files prefixed by "y_" are corresponding labels.
%
%neibour: the number of neighbours on both sides of the original one.
%       This is our inner testing parameter. 
%       In this paper, we set 'neighbour=0', and don't change it. 
%cell: name of the cell line/tissue to be loaded. 


%output:
%The output data is saved in the default directory 
%'Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200'
%and named 'XXX_neibour0.mat', where 'XXX' is the name of the cell/tissue. The matlab file
%'XXX_neibour0.mat' contains two variables, namely 'seqs' and 'states'
%
%seqs:          Features. 'seqs' is a cell vector containing all features of training set
%               , and each cell is a matrix indicating the features of an ehnhancer or a non-enhancer.
%               The length of the cell vector equals the number of ehnhancers and 
%               non-enhanchers in the training set.  For example, if
%               "chr1 1000 2200" is enhancer, then in its corresponding cell
%               is a D*6 matrix. D is the dimension of the feature(for
%               exmaple 1114).Six is the lenght of "chr1 1000 2200", since
%               (2200-1000)/200=6. 200 is the resolution,which means the genome is divided into 200 bp intervals.
%states:        Lables. Similar with 'seqs'. 'states' is a cell vector containing all lables of training set
%               , and each cell is a vector indicating the lables of an ehnhancer or a non-enhancer (1 for enhancer, 2 for non-enhancer).
%               The length of the cell vector equals the number of ehnhancers and 
%               non-enhanchers in the training set. For example, if
%               "chr1 1000 2200" is enhancer, then its corresponding cell
%               is a [1 1 1 1 1 1] vector.



if nargin<1
	neibour=0;      %the number of neighbours on both sides of the original one.
                    %This is our inner testing parameter. 
                    %In this paper, we set 'neighbour=0', and don't change it. 
end
if nargin<2
	cell='H1_cell_line';
end

disp(num2str(neibour));
disp(cell);

num_of_seq=20000;    %maximum number of sub-sequences on one chromosome.

dir_in=strcat('Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200/',cell);    
seqs=[];
states=[];

for i=1:23   % chromosome
	chr=strcat('chr',num2str(i));
	if i==23
		chr='chrX';
    end
    
    disp(chr);
    
	for j=1:num_of_seq
		if ~exist(strcat(dir_in,'/','x_',chr,'_seq',num2str(j),'.matrix'))
			break;
		end
		
		x=load( strcat(dir_in,'/','x_',chr,'_seq',num2str(j),'.matrix') );
		y=load( strcat(dir_in,'/','y_',chr,'_seq',num2str(j),'.matrix') );
		
		if size(x,1)~=size(y,1);
			error('number of sample not match!');
		end
		
		train_x = format_x_locally(x,neibour);
		train_y = y;
		seqs{end+1} = train_x';
		states{end+1} = train_y';
	end
end


save(strcat(dir_in,'/../',cell,'_neibour',num2str(neibour),'.mat'),'seqs','states','-v7.3');


end

