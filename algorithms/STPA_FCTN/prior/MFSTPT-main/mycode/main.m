%this matlab code implements the MFSTPT model for infrared small target 
% detection.
%
% Reference:
% Hu Y, Ma Y, Pan Z, et al. Infrared Dim and Small Target Detection from 
% Complex Scenes via Multi-Frame Spatial–Temporal Patch-Tensor Model[J]. 
% Remote Sensing, 2022, 14(9): 2234.
%
%
% Written by Yapeng Ma 
% 2022-5-05
clear
clc
addpath('functions/')
datafile_path='.\img\';%文件
dir_res='.\result\';%结果
num=20;
%% 模型参数
patchSize = 60;
slideStep = 60;
lambdaL =1;  %tuning   
opts.mu=1*1e-3;
opts.epsilon=1;
opts.rho=1.1;
%% 
fprintf('存储位置: %s\n', dir_res);
img_path_list = dir(fullfile(datafile_path,'*.bmp'));%获取该文件夹中所有bmp格式的图像  
files = length(img_path_list);%获取图像总数量 
files_name =sort_nat({img_path_list.name});%重新排序
 for i = 1:num+4
    newname=fullfile(datafile_path, files_name{i});
    img = imread(newname);
    nn = ndims(img);
    if nn==3
       img= rgb2gray(img);
    end
    D(:,:,i) = double(img);
    end
    for i =3:num+2%+2是保证最后一张图片有时序信息
    %% construct spatial-temporal tensor
        Dpart = D(:,:,i-2:i+2);
        pad_D=padd(Dpart,slideStep);   
        for ii = 1:5
            img = pad_D(:,:,ii);
            [lambda1, lambda2] = structure_tensor_lambda(img, 3);
             %      step 2: calculate corner strength function
            cornerStrength = (((lambda1.*lambda2)./(lambda1 + lambda2)));
            %      step 3: obtain final weight map
             EValue = (lambda1-lambda2);
             priorWeight=mat2gray(real((real(cornerStrength.^1.5)).*(EValue).^(1/2.5)));
             Wetensor(:,:,ii)=priorWeight;
        end
        tenD = gen_patch(pad_D, patchSize, slideStep);
        tenW = gen_patch(Wetensor, patchSize, slideStep);        
%% low-rank and sparse recovery
    fprintf('%s\n',datafile_path);
    fprintf('%d/%d: %s\n', files, i, files_name{i});
    Nway = size(tenD);
    tenB1=zeros(patchSize,patchSize,Nway(4));
    tenT1=zeros(patchSize,patchSize,Nway(4));
    for k=1:Nway(4)
        tenD1=tenD(:,:,:,k);
        tenW1=tenW(:,:,:,k);
        if mod(k,15)==0
            fprintf('第%d/%d块\n',k,Nway(4));
        end
        lambda = lambdaL / sqrt(max(Nway(1),Nway(2))*Nway(3)); 
        [tenB,tenT] = MFSTPT(tenD1,lambda,tenW1,opts);
        tenB1(:,:,k)=tenB(:,:,23);%每一块作为一纵向排列输入进去，得出的结果中中间的块为要求的
        tenT1(:,:,k)=tenT(:,:,23);
    end
        tarImg = res_patch(tenT1, Dpart(:,:,1), patchSize, slideStep);
        %backImg = res_patch_ten_mean(tenB1, D11(:,:,1), patchSize, slideStep);
        %sum(sum(tarImg))
        maxv = max(max(double(img)));
        Target = uint8( mat2gray(tarImg)*maxv );
        imshow(Target,[]),title('Target image')
        %subplot(133),imshow(A,[]),title('Background image') 
        imwrite(Target, fullfile(dir_res,[files_name{i}]));
        %imwrite(A, fullfile(B1,[num2str(i)  '.bmp']));
    end