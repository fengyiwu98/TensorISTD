%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Guanghui Wang
% Email: wangguang147@qq.com
% Copyright:  University of Electronic Science and Technology of China
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

clear;clc;close all;
addpath 'function'
%%%%%%%��ȡͼ���·��%%%%%%%%
% ͼ�������ļ���·��
file_path = 'D:\�ж���\code\benchmark�ϼ�\benchmark_release(wyl)\dataset\data\data13\';
% �W??��ȡͼƬ�Ğ@?�㣨�������ҪС�ڵ����ļ�����ͼƬ�@?�㣩
img_num = 396;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
png_res_base_path='D:\�ж���\code\benchmark�ϼ�\benchmark_release(wyl)\png_results\data13\NPSTT\';
%% ���������Z
%%%%%%%%%���²���������%%%%%%%%%%
patchSize = 70;  % patch �Ĵ���
slideStep = 70;  % �������ڵĲ���
frame = 3;       % ǰ���frame��
patchNum = 3;    % patchNum * patchNum دpatches�߸�patchد?����
curframe = 20;   % ��ǰ֡λ�Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tem = fix(patchNum/2); % ����padding����ʱ����������ȡ����

%% 1. ��ȡ��??ͼƬ��patch�}
% ��ȡ�任�ĲΔ�
 path1 = [file_path,num2str(1,'%04d'),  '.bmp']; 
%path1 = [file_path,num2str(1),'.png'];
img1 = imread(path1); 

if ndims( img1 ) == 3
    img1 = rgb2gray( img1 );
end
tic
%����ͬһ�����еĴ�Сد?�� ������Щ����
[imgHei, imgWid] = size(img1);
% ��???��
rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;  % patch����
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;  % patch����
% ��patch�����Ͻ�����
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];
% ������ų�??ͼƬ��Patch��cell
% extend_curcel Ϊ�˽���߽�����
curcel = cell(rowPatchNum,colPatchNum,img_num);  % ����*����*֡��
extend_curcel = cell(rowPatchNum+2*tem,colPatchNum+2*tem,img_num);  % padding

for i = 1:img_num
    path = [file_path,num2str(i,'%04d'),  '.bmp'];

   % path = [file_path,num2str(i),'.png'];
    img = imread(path); 
   
    %%��ʽ����
    if ndims( img ) == 3
        img = rgb2gray( img );
    end
    img = double(img);
       
    %% ������ͼ��ĳ�??patch���뵽cellد 
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

%% ����NPSTT
tensorNum = patchNum*patchNum*(frame*2+1);  % 7֡ͼ��ÿ֡����patch ��?dim(3)=63
% ��ǰ����63�����е�λ��
% curpatch = patchNum*patchNum*2*frame+patchNum*tem+tem+1;  %���ˣ�����
curpatch = patchNum*patchNum*frame+patchNum*tem+tem+1; 
% Ŀ���뱳�������Ĵ�С���뵱ǰ֡����patch�ѵ��ɵ�������Сد??
curframetartensor = zeros(patchSize, patchSize, rowPatchNum*colPatchNum);
curframebartensor = zeros(patchSize, patchSize, rowPatchNum*colPatchNum);
t = 0;

for row = tem+1:rowPatchNum+tem  % ѡ����ǰ�飨����??padding���֣�
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

        % lambdaΪƽ�����ϡ������Ĳ���
        lambda = 1/sqrt(patchNum*patchNum*(1+2*frame)*patchSize);
        [tenBack,tenTar] = tCSVT(tmp_curpatchtensor,lambda);
     
        t = t+1;       
        curframetartensor(:,:,t) = tenTar(:,:,curpatch);
        curframebartensor(:,:,t) = tenBack(:,:,curpatch);           
    end
end

 path2 = [file_path,...
         num2str(curframe,'%04d'),  '.bmp']; 
%path2 = [file_path,num2str(curframe),'.png'];
curImg = imread(path2);   

%%��ʽ����
if ndims( curImg ) == 3
    curImg = rgb2gray( curImg );
end
curImg = double(curImg);

tarImg = res_patch_ten(curframetartensor, curImg, patchSize, slideStep);
backImg = res_patch_ten(curframebartensor, curImg, patchSize, slideStep);
toc
figure,
subplot(131),imshow(curImg,[]),title('Ori image');
subplot(132),imshow(tarImg,[]),title('Target image');
subplot(133),imshow(backImg,[]),title('Back image');
