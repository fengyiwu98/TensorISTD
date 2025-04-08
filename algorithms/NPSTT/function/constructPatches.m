function [extend_curcel,rowPatchNum,colPatchNum] = constructPatches(file_path,img_num,patchSize,slideStep,patchNum)

% 将序列图像中的每一帧划分为patches储存在一个cell中

% file_path - 原序列存储文件夹
% img_num - 原序列帧数
% patchSize - patch边长
% slideStep - 滑动步长
% patchNum - 邻域由 patchNum * patchNum 个patch构成

% extend_curcel - 将每帧图像划分为 rowPatchNum * colPatchNum 个patch，
%                 然后上下左右各padding tem 个patch，用cell保存

tem = fix(patchNum/2); % padding的patch数

%% 随便找一张图，获取变换的参数
% 同一个序列的大小一样,所以能共享这些参数
path1 = [file_path, num2str(1),  '.png'];
img1 = imread(path1); 

if ndims( img1 ) == 3
    img1 = rgb2gray( img1 );
end

[imgHei, imgWid] = size(img1);
% 行、列数
rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;
% 各patch的左上角坐标
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];
% 创建存放所有帧图片patch的cell
curcel = cell(rowPatchNum,colPatchNum,img_num);  % 行数*列数*总帧数
extend_curcel = cell(rowPatchNum+2*tem,colPatchNum+2*tem,img_num);  % padding


%% 把序列图像的所有patch放入到cell中
for i = 1:img_num  % 遍历每一帧
    
    path = [file_path, num2str(i),  '.png'];  
    img = imread(path); 
   
    if ndims( img ) == 3
        img = rgb2gray( img );
    end
    img = double(img);
    
    % 把第i帧的所有patch放入extend_curcel{xx,xx,i}中
    for row = 1:rowPatchNum
        for col = 1:colPatchNum
            % 在第i张图像中截取第row行第col列的patch
            tmp_patch = img(rowPosArr(1,row) : rowPosArr(1,row) + patchSize - 1, colPosArr(1,col) : colPosArr(1,col) + patchSize - 1);
            curcel{row,col,i} = tmp_patch;
            extend_curcel{row+tem,col+tem,i} = tmp_patch;
        end
    end
    
    % 扩展padding（复制相邻像素）
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