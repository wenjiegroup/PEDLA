function [ record ] = B2o1_merge_performance_of_training_or_test_set( PI,A,B_DNN,test_or_training )
%   Detailed explanation goes here

%   This function evaluates the performance of a trained PEDLA  
%   on 22 training cells/tissues or 20 independent test cells/tissues.
%   The order of the 22 training cells/tissues and 20 independent test cells/tissues
%   can be seen below, which is consistent with the order in the supplementary table of our paper.


%input:
%PI:    initial state probability vecort of HMM of a trained PEDLA 
%A: state transition probability matrix of HMM of a trained PEDLA 
%B_DNN: all model parameters of DNN of a trained PEDLA 
%test_or_training: should be either 'test' or 'training' which indicates evaluation on 20 independent test cells/tissues or 22 training cells/tissues.
%                  defualt: 'test'

%output:
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



if nargin<4
	test_or_training='test';    % test or training
end
if ~strcmp(test_or_training,'test') && ~strcmp(test_or_training,'training')
    error('test_or_training must be test or training!');
end

if strcmp(test_or_training,'test')
    cell_lines={'Hepg2','Dnd41','Nhek','Nhlf','Hsmmt','spleen','gastric','adrenal_gland','small_intestine' ...
        ,'CD4_primary_cells','HUES64_cell_line','heart_left_ventricle','pancreatic_islets','iPS_DF_19_11_cell_line' ...
        ,'neurosphere_cultured_cells_ganglionic_eminence_derived','heart_aorta','penis_foreskin_keratinocyte_primary_cells' ...
        ,'peripheral_blood_mononuclear_primary_cells','hESC-derived_CD56_ectoderm_cultured_cells','brain_hippocampus_middle'};   %20 independent test cells/tissues
elseif strcmp(test_or_training,'training')
    cell_lines={'H1_cell_line','Gm12878','K562','Helas3','Huvec','Monocd14ro1746','Hmec','Hsmm','Nha','Nhdfad','Osteobl','IMR90_cell_line','liver','pancreas','lung','heart_right_atrium','esophagus','psoas_muscle','ovary','sigmoid_colon','thymus','penis_foreskin_fibroblast_primary_cells'};  %22 training cells/tissues
else
    error('test_or_training must be test or training!');
end


record=[];
for kkk=1:length(cell_lines)
    
    [ performance_test ] = B2o2_performance_of_learnt_model_for_one_cell( 0,PI,A,B_DNN,cell_lines{kkk} );
    record=[record performance_test];
end


end

