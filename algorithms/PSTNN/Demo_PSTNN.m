% This matlab code implements the infrared small target detection model
% based on partial sum of the tensor nuclear norm.
% 
% Reference:
% Zhang, L.; Peng, Z. Infrared Small Target Detection Based on Partial Sum 
% of the Tensor Nuclear Norm. Remote Sens. 2019, 11, 382.
%
% Written by Landan Zhang 
% 2019-2-24
clc;
clear;
close all;

%addpath('functions/')
%addpath('tools/')
%saveDir = 'results/';
%% imgpath = 'images/';
%imgpath = 'E:\3_Code\1_Small_target_detection\1_mytest_IRdetection\data\sirst\image\';
%imgDir = dir([imgpath '*.png']);
addpath('./functions/')
addpath('./tools/')

datas = {
    %'sequence1','sequence4','sequence7','sequence9','sequence11','sequence12'
    %'noise3','easy3','noises3','noise_t3','noise_c3','_noise3'
    %'noise5','easy5','noises5','noise_t5','noise_c5','_noise5'
    %'noise6','easy6','noises6','noise_t6','noise_c6','_noise6'
    'noise_t5'
    }
for i = 1:length(datas)
data = datas{i};
%saveDir = 'D:\Tensor\benchmark_release(lsm)\benchmark_release(lsm)\all_result\\sequence7\PSTNN_1e3\';
%imgpath = 'D:\Tensor\benchmark_release(lsm)\benchmark_release(lsm)\dataset\data\sequence7\';
imgpath = ['D:\Tensor\benchmark_release(lsm)\benchmark_release(lsm)\dataset\data\',data,'\'];
saveDir = ['D:\Tensor\benchmark_release(lsm)\benchmark_release(lsm)\all_result\',data,'/PSTNN_1e4/'];
imgDir = dir([imgpath '*.bmp']);
time_path = 'D:\Tensor\benchmark_release(lsm)\benchmark_release(lsm)\time_results\';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% patch parameters
patchSize = 40;
slideStep = 40;
lambdaL = 0.7;  %tuning
total_time = 0;  % 确保 total_time 在循环之前初始化

% flag 是否显示结果
flag = 1;
time = 0;%计算程序耗时

len = length(imgDir);
for i=1:len
%     img = imread([imgpath imgDir(i).name]);
    img = imread([imgpath sprintf('%04d.bmp', i)]);
    %if flag
     % figure,subplot(131)
      %imshow(img),title('Original image')
    %end
    if ndims( img ) == 3
        img = rgb2gray( img );
    end
    img = double(img);
    
    %tic
    %% constrcut patch tensor of original image
    tenD = gen_patch_ten(img, patchSize, slideStep);
    [n1,n2,n3] = size(tenD);  
    
    %% calculate prior weight map
    %      step 1: calculate two eigenvalues from structure tensor
    [lambda1, lambda2] = structure_tensor_lambda(img, 3);
    %      step 2: calculate corner strength function
    cornerStrength = (((lambda1.*lambda2)./(lambda1 + lambda2)));
    %      step 3: obtain final weight map
    maxValue = (max(lambda1,lambda2));
    priorWeight = mat2gray(cornerStrength .* maxValue);
    %      step 4: constrcut patch tensor of weight map
    tenW = gen_patch_ten(priorWeight, patchSize, slideStep);
    
    %% The proposed model
    lambda = lambdaL / sqrt(max(n1,n2)*n3); 
    [tenB,tenT] = trpca_pstnn(tenD,lambda,tenW); 
    
    %% recover the target and background image
    tarImg = res_patch_ten_mean(tenT, img, patchSize, slideStep);
    backImg = res_patch_ten_mean(tenB, img, patchSize, slideStep);

    tic
    [tenB,tenT] = trpca_pstnn(tenD,lambda,tenW);    
    toc
    time=toc+time
   % t = toc
   % time = time + t;
    maxv = max(max(double(img)));
    tarImg = tarImg/max(tarImg(:));
    backImg = backImg/max(backImg(:));
    E = uint8( mat2gray(tarImg)*maxv );
    A = uint8( mat2gray(backImg)*maxv );
    
    %if flag
        %subplot(132),imshow(E,[]),title('Target image')
        %subplot(133),imshow(A,[]),title('Background image')
    %end
    % save the results
%     imwrite(E, [saveDir 'target/' imgDir(i).name]);
%     imwrite(A, [saveDir 'background/' imgDir(i).name]);
      imwrite(E, [saveDir  sprintf('%04d.png', i)]);
      %imwrite(A, [saveDir 'background\PSTNN_background_' num2str(i) '.png']);
   % elapsed_time = toc;  % 获取当前图片的处理时间
    %total_time = total_time + elapsed_time;  % 累加到总时间
    
    avg_time = time / 100;  % 平均每次处理时间
    txt_path = [time_path, 'PSTNN_1e4_time.txt'];
    fid = fopen(txt_path, 'a');
    time_name = datestr(now, 'yy-mm-dd_HH-MM-SS');
    %fprintf(fid, '%s\n', [time_name, ':']);
    fprintf(fid, '%s\n', ['PSTNN in ',data,' took ', num2str(avg_time, '%.4f'), ' seconds.']);
    fclose(fid);
end
end  