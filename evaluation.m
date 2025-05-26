%% Benchmark of Infraed Samll Target Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author:
% Email: wufengyi98@163.com
% Copyright:  University of Electronic Science and Technology of China
% Date: 2025/04/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.
fclose('all'); close all; clear all; clc;

algo_base_path = '.\algorithm\';
data_base_path = '.\dataset/data\';
res_base_path =  '.\all_results\';
time_path = '.\time_results\';
mat_base_path = '.\mat_results\';
figure_base_path = '.\fig_results\';
txt_base_path =  './dataset/anno/';
mask_base_path = './';
addpath('./utils/');

%% user configs
eval_algo_names = ...
    {
     'STPA_FCTN'
    };

eval_data_names = ...
    {
     'sequence1'
    };

img_types = {'*.jpg', '*.bmp', '*.png'};
preimg_type = '*.png';
thres_num = 100;
radius = 1; % radius of object [7x7]
x_axis_ratio = 1e-5;
FPR_thres = 1;

%% step 1: Use all the evaluation algorithms to get the result plots for all the evaluation datasets present in . /result in the mat file
get_algo_result(eval_algo_names, eval_data_names, ...
     img_types, algo_base_path, data_base_path, res_base_path, time_path );
%Multiframe algorithms(STT, NPSTT)
%get_algo_result_multiframe(eval_algo_names, eval_data_names, ...
     %img_types, algo_base_path, data_base_path, res_base_path, time_path );

%% step 2: Combining the result plots from step 1 and . GT under /dataset/ann/ to get the roc sequence results
get_curves(eval_algo_names, eval_data_names, thres_num, radius, res_base_path, ...
     mat_base_path, txt_base_path, mask_base_path, preimg_type);

%% step 3: Calculating metrics and plotting 3DROC
curves_drawer(1, eval_algo_names, eval_data_names, figure_base_path, mat_base_path, x_axis_ratio, FPR_thres);

%% step 4: Calculate SCRG gain,CG, BSF, BSR
measure_calculator(eval_algo_names, eval_data_names, data_base_path, res_base_path, ...
    mat_base_path, txt_base_path, img_types, preimg_type);
