%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Guanghui Wang
% Email: wangguang147@qq.com
% Copyright:  University of Electronic Science and Technology of China
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

clear all;clc;close all;
addpath 'function'
addpath 'lib'
%%%%%%%��ȡͼ���·��%%%%%%%%
% ͼ�������ļ���
exter='.bmp';


%% ���������Z
%%%%%%%%%���²������޸�%%%%%%%%%%
patchSize = 70;  % patch 
slideStep = 70;  % �������ڵĲ���
frame = 3;       % ǰ���frame
patchNum = 3; 
tem = fix(patchNum/2); % ����padding����ʱ����������ȡ����
Length=50;
for tt=3:3
seq=['_noise' num2str(tt) '\'];
img_path = ['D:\dataset\data\',seq];
des_path=['D:\all_result\',seq,'\NPSTT\'];
if ~exist(des_path, 'dir')
    mkdir(des_path);
end


img_num = Length;
initframe=4;
endframe=img_num-3;
 path1 = [img_path,num2str(1,'%04d'),  exter]; 
%path1 = [file_path,num2str(1),'.png'];
img1 = imread(path1); 

if ndims( img1 ) == 3
    img1 = rgb2gray( img1 );
end

%����ͬһ�����еĴ�С ������Щ����
[imgHei, imgWid] = size(img1);
% ����
rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;  % patch����
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;  % patch����
% ��patch�����Ͻ�����
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];
% �������ͼƬ��Patch��cell
% extend_curcel Ϊ�˽���߽�����
curcel = cell(rowPatchNum,colPatchNum,img_num);  % ����*����*֡��
extend_curcel = cell(rowPatchNum+2*tem,colPatchNum+2*tem,img_num);  % padding

for i = 1:img_num
    path = [img_path,num2str(i,'%04d'),  exter];

   % path = [file_path,num2str(i),'.png'];
    img = imread(path); 
   
    %%��ʽ����
    if ndims( img ) == 3
        img = rgb2gray( img );
    end
    img = double(img);
tic
       
    %% ������ͼ���patch���뵽cellد 
    for row = 1:rowPatchNum
        for col = 1:colPatchNum
            % �ڵ�i��ͼ���н�ȡ��row�е�col�е�patch
            tmp_patch = img(rowPosArr(1,row) : rowPosArr(1,row) + patchSize - 1, colPosArr(1,col) : colPosArr(1,col) + patchSize - 1);
            curcel{row,col,i} = tmp_patch;
            extend_curcel{row+tem,col+tem,i} = tmp_patch;
        end
    end
    
    % ��չ�������������أ�
    % �Ľ�
    for t1 = 1:tem
        for t2 = 1:tem
            extend_curcel{t1,t2,i} = curcel{1,1,i};
        end
    end
    for t1 = 1:tem
        for t2 = tem+colPatchNum+1:tem+colPatchNum+tem
            extend_curcel{t1,t2,i} = curcel{1,colPatchNum,i};
        end
    end
    for t1 = tem+rowPatchNum+1:tem+rowPatchNum+tem
        for t2 = 1:tem
            extend_curcel{t1,t2,i} = curcel{rowPatchNum,1,i};
        end
    end
    for t1 = tem+rowPatchNum+1:tem+rowPatchNum+tem
        for t2 = tem+colPatchNum+1:tem+colPatchNum+tem
            extend_curcel{t1,t2,i} = curcel{rowPatchNum,colPatchNum,i};
        end
    end
    % �ı�
    for t1 = 1:tem
        for t2 = tem+1:tem+colPatchNum
            extend_curcel{t1,t2,i} = curcel{1,t2-tem,i};
        end
    end
    for t1 = tem+rowPatchNum+1:tem+rowPatchNum+tem
        for t2 = tem+1:tem+colPatchNum
            extend_curcel{t1,t2,i} = curcel{rowPatchNum,t2-tem,i};
        end
    end
    for t1 = tem+1:tem+rowPatchNum
        for t2 = 1:tem
            extend_curcel{t1,t2,i} = curcel{t1-tem,1,i};
        end
    end
    for t1 = tem+1:tem+rowPatchNum
        for t2 = tem+colPatchNum+1:tem+colPatchNum+tem
            extend_curcel{t1,t2,i} = curcel{t1-tem,colPatchNum,i};
        end
    end
end

for curframindex=initframe:endframe
   curframe =  curframindex;
%% ����NPSTT
tensorNum = patchNum*patchNum*(frame*2+1);  % 7֡ͼ��
% ��ǰ����63�����е�λ��
% curpatch = patchNum*patchNum*2*frame+patchNum*tem+tem+1;  
% curpatch = patchNum*patchNum*frame+patchNum*tem+tem+1; 
% Ŀ���뱳�������Ĵ�С���뵱ǰ֡����patch�ѵ��ɵ�������С
curframetartensor = zeros(patchSize, patchSize, rowPatchNum*colPatchNum);
curframebartensor = zeros(patchSize, patchSize, rowPatchNum*colPatchNum);
t = 0;

for row = tem+1:rowPatchNum+tem  % ѡ����ǰ��
    for col = tem+1:colPatchNum+tem
        
        % ��ȡ��Ӧ��ǰ��֡�ĵ�ǰ���tensor
        k = 0;
        tmp_curpatchtensor = zeros(patchSize, patchSize, tensorNum);        
        for frameNum = curframe-frame:curframe+frame  % �뵱ǰ֡��ص�ǰ��֡
            for patchRow = -tem:tem  % ��ǰ�鼰��ռ���Χ��
                for patchCol = -tem:tem
                    k = k+1;
                    tmp_curpatchtensor(:,:,k) = extend_curcel{row+patchRow,col+patchCol,frameNum};
                end
            end           
        end
        
        lambda = 0.5/sqrt(patchNum*patchNum*(1+2*frame)*patchSize);
         [tenBack,tenTar] = tCSVT(tmp_curpatchtensor,lambda);
        t = t+1;       
        curframetartensor(:,:,t) = tenTar(:,:,curpatch);
        curframebartensor(:,:,t) = tenBack(:,:,curpatch);
      if curframindex==initframe 
            for index=1:3
                L(:,:,t,index)=  tenTar(:,:,(index-1)*patchNum*patchNum+patchNum*tem+tem+1);
            end
        end
         if curframindex==endframe 
            for index=1:3
                L(:,:,t,index)=  tenTar(:,:,(index+4-1)*patchNum*patchNum+patchNum*tem+tem+1);
            end
        end
    end
end

 path2 = [img_path,...
         num2str(curframe,'%04d'),  exter]; 
%path2 = [file_path,num2str(curframe),'.png'];
curImg = imread(path2);   

%%��ʽ����
if ndims( curImg ) == 3
    curImg = rgb2gray( curImg );
end
curImg = double(curImg);

tarImg = res_patch_ten(curframetartensor, curImg, patchSize, slideStep);
%backImg = res_patch_ten(curframebartensor, curImg, patchSize, slideStep);
res_path = [des_path, num2str(curframe, '%04d'), '.png'];
imwrite(tarImg / (max(tarImg(:)) + eps), res_path);

        if curframindex==initframe 
            for index=1:3
                tarImg = res_patch_ten(L(:,:,:,index), curImg, patchSize, slideStep);
                res_path = [des_path,num2str(index, '%04d'), '.png'];
                imwrite(tarImg / (max(tarImg(:)) + eps), res_path);
            end
        end
        if curframindex==endframe 
            for index=1:3
                tarImg = res_patch_ten(L(:,:,:,index), curImg, patchSize, slideStep);
                res_path = [des_path, num2str(index+endframe, '%04d'), '.png'];
                imwrite(tarImg / (max(tarImg(:)) + eps), res_path);
            end
        end
 
end
toc
%time_name = datestr(now, 'yy-mm-dd_HH-MM-SS');
% res_base_path='D:\code\tensor_code\shiyan\benchmark_release(sea))_2\time_test_results\';
% txt_path = [res_base_path, 'time_NPSTT_sequence',  num2str(tt) , '.txt'];
% fid = fopen(txt_path, 'a');
% times = 1.0 * toc/100 ;
% disp( num2str(times, '%.4f'));
% fprintf(fid, '%s\n', num2str(times, '%.4f'),'times.');
% fclose(fid);
end

