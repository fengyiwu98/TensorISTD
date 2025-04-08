%% Benchmark of Infraed Samll Target Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author:
% Email: 
% Copyright:  University of Electronic Science and Technology of China
% Date: 2025/04/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.
fclose('all'); close all; clear all; clc;

algo_base_path = '.\algorithm\';
data_base_path = '.\dataset/data\';
res_base_path =  '.\all_result\';
time_path = '.\time_results\';
mat_base_path = '.\mat_results\';
figure_base_path = '.\fig_results\';
txt_base_path =  './dataset/anno/';
mask_base_path = './';
addpath('./utils/');

%% user configs
eval_algo_names = ...
    { 
    %'$a=1$','$a=2$','$a=3$'
    %'$I_{ws}=30$','$I_{ws}=40$','$I_{ws}=50$','$I_{ws}=60$','$I_{ws}=70$','$I_{ws}=80$'
    %'$I_{t}=5$','$I_{t}=10$','$I_{t}=15$','$I_{t}=20$','$I_{t}=25$',
    'STPA-FCTN','WSWTNN-PnP','3D-STPM','RPCANet','MSHNet','5-D-STFC','4-D-TR','SDD','NPSTT','STT','ASTTV-NTLA','PSTNN','IPI','MPCM'
    %'WSWTNN-PnP','RPCANet','MSHNet'
    %'5-D-STFC'
    %'STPA-FCTN','3D-STPM','4-D-TR','SDD','NPSTT','STT','ASTTV-NTLA','PSTNN','IPI','MPCM'
    %'STPA_FCTN_bb_1e4','STPA_FCTN_a2_1e4','STPA_FCTN_a3_1e4'
    };

eval_data_names = ...
    {
     %'sequence1','sequence4','sequence7','sequence9','sequence11','sequence12'
     %'noise3','noise5','noise6','easy3','easy5','easy6','noises3','noises5','noises6','noise_t3','noise_t5','noise_t6','noise_c3','noise_c5','noise_c6','_noise3','_noise5','_noise6'
     'sequence11'
    };

img_types = {'*.jpg', '*.bmp', '*.png'};
preimg_type = '*.jpg';
thres_num = 100;
radius = 1; % radius of object [7x7]
x_axis_ratio = 1e-4;
FPR_thres = 1;

%% step 1: 使用所有评测算法得到对于所有评测数据集的结果图，存在./result下的mat文件
%%get_algo_result(eval_algo_names, eval_data_names, ...
     %%img_types, algo_base_path, data_base_path, res_base_path, time_path );

%% step 2: 结合step 1的结果图和./dataset/ann/下的GT，得到roc序列结果
%get_curves(eval_algo_names, eval_data_names, thres_num, radius, res_base_path, ...
     %mat_base_path, txt_base_path, mask_base_path, preimg_type);


%% step 3: 结合step 2中的序列结果，绘制出结果图，图数量和评测数据集数量一致
curves_drawer(1, eval_algo_names, eval_data_names, figure_base_path, mat_base_path, x_axis_ratio, FPR_thres);

%% step 4: 计算SCRG增益和BSF
measure_calculator(eval_algo_names, eval_data_names, data_base_path, res_base_path, ...
    mat_base_path, txt_base_path, img_types, preimg_type);
