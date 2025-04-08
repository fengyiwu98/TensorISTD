clc;
clear;
close all;

addpath('./functions/')
%img_path = 'D:\code\tensor_code\shiyan\benchmark_release(sea))_2\dataset\data\easy6\';
%des_path='./results/Seq1/';
%des_path='D:\code\tensor_code\shiyan\comparision\FCTN_prior_2\all_result\easy6\';
patchSize = 70;
slideStep = 70;
start_idx = 1;
length = 100;
L = 15;
for seq=1:1
 time=0
 
data=[num2str(seq) '\'];
img_path=['D:\code\tensor_code\shiyan\benchmark_release(sea))_2\dataset\data\sequence',data];
des_path=['D:\code\tensor_code\shiyan\benchmark_release(sea))_2\time_test_results\sequence',data,'\4_D_TR\'];
%prior_result_path=['D:\code\tensor_code\shiyan\comparision\FCTN_prior_2\prior_result\STLDM\sequence',data];
if ~exist(des_path, 'dir')
    mkdir(des_path);
end
for i = 1:floor(length/L)
    for j=1:L
        cur_idx = (i-1) * L + j + start_idx - 1;
        img = imread( [img_path num2str(cur_idx,'%04d') '.bmp'] );
        if ndims(img) == 3
            img = rgb2gray(img);
        end
        img = im2double(img);
         D = gen_patch_ten(img, patchSize, slideStep);  
         tensor_4D(:,:,j,:)=D;
    end
    
    tic
    [tenBack,tenTar] = TR_RPCA(tensor_4D);    
    toc
    time=toc+time
     for j = 1: L
        cur_idx = (i-1) * L + j + start_idx - 1;
        
        tarImg = res_patch_ten_mean (tenTar(:,:,j,:), img, patchSize, slideStep);
        bacImg = res_patch_ten_mean (tenBack(:,:,j,:), img, patchSize, slideStep);
        
        T = tarImg / (max(tarImg(:)) + eps);
        B = bacImg / (max(bacImg(:)) + eps);
        
        %tar_res_path = [des_path,'TR_target\', num2str(cur_idx,'%04d'), '.png'];
        tar_res_path = [des_path,num2str(cur_idx,'%04d'), '.png'];
        back_res_path = [des_path,'TR_background\', num2str(cur_idx,'%04d'), '.png'];
        
        imwrite(T, tar_res_path);
        %imwrite(B, back_res_path);
    end   
end
    res_base_path='D:\code\tensor_code\shiyan\benchmark_release(sea))_2\time_test_results\';
    txt_path = [res_base_path, 'time_4D-TR_sequence',  num2str(seq) , '.txt'];
    fid = fopen(txt_path, 'a');
    times = 1.0 * time/100 ;
    disp( num2str(times, '%.4f'));
    fprintf(fid, '%s\n', num2str(times, '%.4f'),'times.');
    fclose(fid);

end