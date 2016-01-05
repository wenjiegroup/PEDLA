function [ ] = P2o_calculate_quantile_of_twentytwo_training_cells( neibour )
%   Detailed explanation goes here

%This function calculates the quantile of all features for twentytwo
%training cells/tissues. In our paper, we use the 100-th quantile(namely
%the maximum) to normlize all features.

%output:
%The output data is saved in the default directory 
%'Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200'
%and named 'q100_of_twentytwo_training_cells_neibour0.mat'


if nargin<1
	neibour=0;      %the number of neighbours on both sides of the original one.
                    %This is our inner testing parameter. 
                    %In this paper, we set 'neighbour=0', and don't change it. 
end

cell_lines={'H1_cell_line','Gm12878','K562','Helas3','Huvec','Monocd14ro1746','Hmec','Hsmm','Nha','Nhdfad','Osteobl','IMR90_cell_line','liver','pancreas','lung','heart_right_atrium','esophagus','psoas_muscle','ovary','sigmoid_colon','thymus','penis_foreskin_fibroblast_primary_cells'};  %twentytwo

seqs_temp=[];
states_temp=[];

for kkk=1:length(cell_lines)
    disp(cell_lines{kkk});
    dir_in=strcat('Data_for_training_and_prediction/Data_for_training_and_prediction-ratio9step200win200/');  
    load(strcat(dir_in,'/',cell_lines{kkk},'_neibour',num2str(neibour),'.mat'),'seqs','states');
    seqs_temp=[seqs_temp seqs];
    states_temp=[states_temp states];
end

dir_out=strcat(dir_in);   % twentytwo_training_cells
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end
seqs=seqs_temp;
states=states_temp;
clear seqs_temp;
clear states_temp;

my_quantile=100;  %the quantile used to scale data
temp=[];
for k = 1:numel(seqs)   
    temp=[temp,seqs{k}];
end
q=quantile(temp',[my_quantile/100])';
temp=[];
save(strcat(dir_out,'/q',num2str(my_quantile),'_of_twentytwo_training_cells_neibour',num2str(neibour),'.mat'),'q','-v7.3');


end

