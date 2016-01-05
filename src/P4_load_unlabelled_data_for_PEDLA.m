function [  ] = P4_load_unlabelled_data_for_PEDLA( neibour )
%   Detailed explanation goes here

%   This function just loads the extracted unlabeled features into a
%   matlab format.


%input:
%The features extracted in a previous step are loaded from sub-directories of the default directory 
%'Data_for_training_and_prediction/unlabelled_data-step200win200/'.
%The files prefixed by "x_" are the features.
%
%neibour: the number of neighbours on both sides of the original one.
%       This is our inner testing parameter. 
%       In this paper, we set 'neighbour=0', and don't change it. 



%output:
%The output data is saved in sub-directories of the default directory 
%'Data_for_training_and_prediction/unlabelled_data-step200win200/'
%and named 'chrN_neibour0.mat', where 'N' is index of the chromosome. The matlab file
%'chrN_neibour0.mat' contains one variable, namely 'unlabelled_seqs'.
%
%unlabelled_seqs:          unlabeled features. 'unlabelled_seqs' is a matrix indicating the features of a whole chromosome.
%               For example, the number of base pairs of chr21 is 48129895, then the matrix
%               is a D*240650 matrix. D is the dimension of the feature(for
%               exmaple 1114). 240650 is the lenght of chr21 in 200-bp resolution, since
%               48129895/200=240650. 200 is the resolution,which means the genome is divided into 200 bp intervals.




if nargin<1
	neibour=0;      %the number of neighbours on both sides of the original one.
                    %This is our inner testing parameter. 
                    %In this paper, we set 'neighbour=0', and don't change it. 
end

dir_in='Data_for_training_and_prediction/unlabelled_data-step200win200';

subdirs=dir(dir_in);

for f=3:numel(subdirs);
    disp(subdirs(f).name);
    
    unlabelled_seqs=[];
    for i=1:23
        chr=strcat('chr',num2str(i));
        if i==23
            chr='chrX';
        end
        disp(chr);
        if exist(strcat(dir_in,'/',subdirs(f).name,'/','x_',chr,'.matrix'),'file');
            x=load( strcat(dir_in,'/',subdirs(f).name,'/','x_',chr,'.matrix') );
        else
            continue;
        end
        

        unlabelled_x=format_x_locally(x,neibour);
        unlabelled_seqs = unlabelled_x';
        save(strcat(dir_in,'/',subdirs(f).name,'/',chr,'_neibour',num2str(neibour),'.mat'),'unlabelled_seqs','-v7.3');

    end  


end

