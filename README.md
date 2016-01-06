------------------------------------------------------------
------------------------------------------------------------
The requirements for running PEDLA are:
MATLAB 8.2(R2013b) 
Perl v5.14.2

The core algorithm of PEDLA is written in MATLAB, and some helpful scripts are written in Perl.

Due to the size limitation of the uploaded file, here we provide only codes of PEDLA and do not provide any data. To run PEDLA without error, first download the PEDLA with codes and example data from the website https://mega.nz/#!rEI0nLCT (these codes are the same as codes here).

------------------------------------------------------------
------------------------------------------------------------
This directory contains all the MATLAB codes and Perl scripts we writed for PEDLA.
The PEDLA consists of two main procedures: 
(a)Training procedure
(b)Prediction procedure

As described in our paper, PEDLA was applied to: (1)train a single-cell model from H1 cell using full 1114-dimensional features, and predict enhancers base on the model; (2)train a multiple-cells/tissues model from 22 training cells/tissues, and predict enhancers in the 22 training cells/tissues and another 20 independent test cells/tissues base on the model.

******************************************************
(1)train a single-cell model from H1 cell using full 1114-dimensional features, and predict enhancers base on the model

(a)Training procedure
The 'A1_training_of_PEDLA_for_H1.m' MATLAB script is the main training procedure. This function learns the model parameters of PEDLA. 
The optimal paramerters we used in our paper are neighbour=0;sizes=[500 500];numepochs=50;cell_line='H1_cell_line'. Since the training procedure will take about several hours and 40 G physical memory, we use neighbour=0;sizes=[10];numepochs=5;cell_line='H1_cell_line' as default input parameters for quick investigation or test. You can omit the input parameter, which will use the default input parameters.
There are detailed instructions and comment lines in the script.

(b)Prediction procedure
The 'A2_prediction_of_PEDLA_for_H1.m' MATLAB script is the main prediction procedure. This function predicts the state (enhancers or non-enhancers) of the unlabeled data based on the learned model of PEDLA in the previous training procedure.
For details, see the instructions and comment lines in the script.

The 'A2o_transform_prediction_of_H1_into_Bed.pl' Perl script tansforms the prediction of PEDLA to the 'bed' format. The fourth column of '.bed' file indicates the state with 'E1' indicating enhancer and 'E2' indicating non-enhancer. 

Examples:
To run PEDLA for a quick test, type "A1_training_of_PEDLA_for_H1()" to train and then type "A2_prediction_of_PEDLA_for_H1()" to predict.
To run PEDLA with optimal paramerters as we used in our paper, type "A1_training_of_PEDLA_for_H1( 0, [500 500], 50, 'H1_cell_line' )" to train and then type "A2_prediction_of_PEDLA_for_H1( 0, [500 500], 50, 'H1_cell_line' )" to predict.

******************************************************
(2)train a multiple-cells/tissues model from 22 training cells/tissues, and predict enhancers in the 22 training cells/tissues and another 20 independent test cells/tissues base on the model.

(a)Training procedure
The 'B1_training_of_PEDLA_for_multiple_cells_and_tissues.m' MATLAB script is the main training procedure. This function learns the model parameters of PEDLA from multiple cells and tissues. The PEDLA for multiple cells and tissues is a general framework of PEDLA and can be used to train on multiple cell lines/tissues. PEDLA training a single-cell model is a special case of PEDLA for multiple cells and tissues where the number of training cells/tissues used for PEDLA is 1.
The optimal paramerters we used in our paper are neighbour=0;sizes=[50];numepochs=150;order_ind=1. Since the training procedure is extremly time-consuming and finishing training on the whole 22 training cells and tissues may take a few days, we use neighbour=0;sizes=[5];numepochs=5;order_ind=1 as default input parameters for quick investigation or test, which can finish training very qucikly. You can omit the input parameter, which will use the default input parameters.
There are detailed instructions and comment lines in the script.

(b)Evaluation procedure
The 'B2_performance_evaluation_of_on_training_or_test_set.m' MATLAB script is the evaluation procedure used in our paper. The evaluation procedure is optional, since the training and prediction procedure of the PEDLA are necessary. This function evaluates the performance of a trained PEDLA model on our 22 training cells/tissues or 20 independent test cells/tissues.
For details, see the instructions and comment lines in the script.
MATLAB scripts 'B2o1_merge_performance_of_training_or_test_set.m' and 'B2o2_performance_of_learnt_model_for_one_cell.m' are sub-functions used by 'B2_performance_evaluation_of_on_training_or_test_set.m'.

(c)Prediction procedure
The 'B3_prediction_of_PEDLA_for_multiple_cells_and_tissues.m' MATLAB script is the main prediction procedure. This function predicts the state (enhancers or non-enhancers) of the unlabeled data based on the learned model of PEDLA in the previous training procedure. This function is almost the same with
 'A2_prediction_of_PEDLA_for_H1.m'.
For details, see the instructions and comment lines in the script.

The 'B3o_transform_prediction_into_Bed.pl' Perl script tansforms the prediction of PEDLA to the 'bed' format. The fourth column of '.bed' file indicates the state with 'E1' indicating enhancer and 'E2' indicating non-enhancer. 

Examples:
To run PEDLA for a quick test, type "B1_training_of_PEDLA_for_multiple_cells_and_tissues()" to train and then type "B3_prediction_of_PEDLA_for_multiple_cells_and_tissues()" to predict.
To run PEDLA with optimal paramerters as we used in our paper, type "B1_training_of_PEDLA_for_multiple_cells_and_tissues( 0, [50], 150, 1 )" to train and then type "B3_prediction_of_PEDLA_for_multiple_cells_and_tissues( 0, [50], 150, 1, 22 )" to predict.



------------------------------------------------------------
------------------------------------------------------------
Preparation of Input Data for PEDLA

The above is our standard training and prediction procedure of PEDLA after the input data are well prepared and transformed into our standard input format of PEDLA. 

******************************************************
Calculation of the RPKM-like values

All the Perl scripts with a name started with 'O' are scripts we used in our paper to calculate the RPKM-like values of features. 

In our paper, we used 9 categories of heterogeneous features, including histone modifications, transcription factors and cofactors, chromatin accessibility, transcription, DNA methylation, CpG island, evolutionary conservation, sequence signatures, and occupancy of TFBS motifs. The calculations of them are almost the same and the corresponding raw data is extremly huge (more than 1 TB), so here we only provide calculation scripts of two categories of features, namely histone modifications and evolutionary conservation.

Scripts 'O1a_tag_count_for_histone_modification.pl', 'O1b_RPKM_calculation_for_histone_modification.pl' and 'O1c_merge_RPKM_of_replicates_for_histone_modification.pl' are the calculations of histone modification features. To save space, here we just provide a example of raw data in the directory 'Raw_Data/Histone_Modification/H1_cell_line'. You can download and put more raw data in the corresponding locations with a same naming mode. Typing the three Perl scripts in order in the command line will run the calculations of histone modification features. 
For details, see the instructions and comment lines in the script.

Scripts 'O2a_tag_count_for_phastCons46way.pl' and 'O2b_RPKM_calculation_for_phastCons46way.pl' are the calculations of evolutionary conservation feature. To save space, here we just provide the raw data of chromosome 21 in the directory 'Raw_Data/phastCons46way'. To run the Perl scripts without error, you need to download the raw data of the other 22 chromosomes (chr1~chr20, chr22 and chrX) to this directory from http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate/. Typing the two Perl scripts in order in the command line will run the calculations of evolutionary conservation feature. 
For details, see the instructions and comment lines in the script.

Script 'O10_combine_data_of_H1.pl' combines 9 categories of heterogeneous data to form matrix as we described in our paper. The total dimensionality of all the heterogeneous data is 1114. 
Script 'O10_combine_data.pl' combines histone modification/phastCons46way data to form matrix as we described in our paper. Histone modification includes H3K4me1|H3K4me3|H3K9me3|H3K27me3|H3K36me3|H3K27ac|input, so the total dimensionality of histone modification/phastCons46way is 8.

Notes:
Because both the raw data and temporary results are extremly huge (much more than 1 TB), here we just provide a example of raw data and don't provide any temporary results.

******************************************************
Transformation into our standard input format of PEDLA

All the MATLAB and Perl scripts with a name started with 'P' are scripts we used in our paper to transform the features and labels (or unlabeled features) into our standard input format of PEDLA after the calculation of the RPKM-like values. 

Runing 'P1_extract_training_feature_and_label_for_PEDLA_of_H1.pl' and 'P2_load_training_feature_and_label_for_PEDLA_of_H1.m' in order will transform the features and labels into our standard input format of PEDLA for training a single-cell model from H1 cell using full 1114-dimensional features, and runing 'P3_extract_unlabelled_data_for_PEDLA_of_H1.pl' and 'P4_load_unlabelled_data_for_PEDLA_of_H1.m' in order will transform the unlabeled features into our standard input format of PEDLA for predicting enhancers base on the model.
For details, see the instructions and comment lines in the script.

Runing 'P1_extract_training_feature_and_label_for_PEDLA.pl' and 'P2_load_training_feature_and_label_for_PEDLA1.m' in order will transform the features and labels into our standard input format of PEDLA for training a multiple-cells/tissues model from 22 training cells/tissues, and runing 'P3_extract_unlabelled_data_for_PEDLA.pl' and 'P4_load_unlabelled_data_for_PEDLA.m' in order will transform the unlabeled features into our standard input format of PEDLA for predicting enhancers in the 22 training cells/tissues and another 20 independent test cells/tissues base on the model.
Runing 'P2o_calculate_quantile_of_twentytwo_training_cells.m' will calculate the quantile of all features for 22 training cells/tissues used to normlize all features.
For details, see the instructions and comment lines in the script.


The standard '.mat' file of PEDLA for training contains two variables, namely 'seqs' and 'states'. 'seqs' are the features. 'seqs' is a cell vector containing all features of training set, and each cell is a matrix indicating the features of an ehnhancer or a non-enhancer. The length of the cell vector equals the number of ehnhancers and non-enhanchers in the training set.  For example, if "chr1 1000 2200" is enhancer, then in its corresponding cell is a D*6 matrix. D is the dimension of the feature(for exmaple 1114). Six is the lenght of "chr1 1000 2200", since (2200-1000)/200=6. 200 is the resolution, which means the genome is divided into 200 bp intervals.
'states' are lables. Similar with 'seqs', 'states' is a cell vector containing all lables of training set, and each cell is a vector indicating the lables of an ehnhancer or a non-enhancer (1 for enhancer, 2 for non-enhancer). The length of the cell vector equals the number of ehnhancers and  non-enhanchers in the training set. For example, if "chr1 1000 2200" is enhancer, then its corresponding cell is a [1 1 1 1 1 1] vector.

The standard '.mat' file of PEDLA for prediction contains one variable, namely 'unlabelled_seqs'.
'unlabelled_seqs' are unlabeled features. 'unlabelled_seqs' is a matrix indicating the features of a whole chromosome. For example, the number of base pairs of chr21 is 48129895, then the matrix is a D*240650 matrix. D is the dimension of the feature(for exmaple 1114). 240650 is the lenght of chr21 in 200-bp resolution, since 48129895/200=240650. 200 is the resolution,which means the genome is divided into 200 bp intervals.





