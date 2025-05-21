clear;clc;close all;

patchSize = 50; 
slideStep = 30; 

frame = 3; 
patchNum = 3; 
tem = fix(patchNum/2);
Length=50;
for tt=6:6
seq=['easy' num2str(tt) '\'];

img_path = ['D:\dataset\data\',seq];
des_path=['D:time_test_results\',seq,'STT\'];
 if ~exist(des_path, 'dir')
    mkdir(des_path);
 end

img_num = Length;
 initframe=4;
 endframe=img_num-3;
   if j>1 && j<6
       path1 = [img_path,...
           num2str(1,'%04d'), '.bmp'];
   elseif j == 9
       path1 = [img_path,...
           num2str(1,'%04d'), '.bmp'];
   else
       path1 = [img_path,...
           num2str(1,'%04d'), '.bmp'];
   end
   
%%
img1 = imread(path1); 

if ndims( img1 ) == 3
    img1 = rgb2gray( img1 );
end

[imgHei, imgWid] = size(img1);

rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;

rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];

curcel = cell(rowPatchNum,colPatchNum,img_num);
extend_curcel = cell(rowPatchNum+2*tem,colPatchNum+2*tem,img_num);

for i = 1:img_num
    if j>1 && j<6
        path = [img_path,...
            num2str(i,'%04d'), '.bmp'];
    elseif j == 9
        path = [img_path,...
            num2str(i,'%04d'), '.bmp'];
    else
        path = [img_path,...
            num2str(i,'%04d'), '.bmp'];
    end
    
    img = imread(path); 
   
    %%
    if ndims( img ) == 3
        img = rgb2gray( img );
    end
    img = double(img);
       
    %% 
    for row = 1:rowPatchNum
        for col = 1:colPatchNum
            tmp_patch = img(rowPosArr(1,row) : rowPosArr(1,row) + patchSize - 1, colPosArr(1,col) : colPosArr(1,col) + patchSize - 1);
            curcel{row,col,i} = tmp_patch;
            extend_curcel{row+tem,col+tem,i} = tmp_patch;
        end
    end
    
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
time=0;
for curframindex = initframe:endframe 
    tic
%% patch tensor
curframe =  curframindex;
tensorNum = patchNum*patchNum*(frame*2+1); 
curpatch = patchNum*patchNum*frame+patchNum*tem+tem+1;
curframetartensor = zeros(patchSize, patchSize, rowPatchNum*colPatchNum);
curframebartensor = zeros(patchSize, patchSize, rowPatchNum*colPatchNum);
t = 0;

for row = tem+1:rowPatchNum+tem  
    for col = tem+1:colPatchNum+tem
        
        k = 0;
        tmp_curpatchtensor = zeros(patchSize, patchSize, tensorNum);        
        for frameNum = curframe-frame:curframe+frame
            for patchRow = -tem:tem
                for patchCol = -tem:tem
                    k = k+1;
                    tmp_curpatchtensor(:,:,k) = extend_curcel{row+patchRow,col+patchCol,frameNum};
                end
            end           
        end
        
        tmp_curpatchtensor1 = tmp_curpatchtensor(:,:,patchNum*patchNum*frame+1:patchNum*patchNum*(frame+1));
        tmp_curpatchtensor2 = tmp_curpatchtensor(:,:,patchNum*patchNum*(frame+1)+1:tensorNum);
        tmp_curpatchtensor(:,:,patchNum*patchNum*frame+1:patchNum*patchNum*2*frame) = tmp_curpatchtensor2;
        tmp_curpatchtensor(:,:,patchNum*patchNum*2*frame+1:tensorNum) = tmp_curpatchtensor1;

        lambda = 1/sqrt(patchNum*patchNum*(1+2*frame*patchSize)+slideStep);
        [tenBack,tenTar] = TRPCA(tmp_curpatchtensor, lambda);  
       
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

    if j>1 && j<6
        path3 = [img_path,...
            num2str(curframe,'%04d'), '.bmp'];
    elseif j == 9
        path3 = [img_path,...
            num2str(curframe,'%04d'), '.bmp'];
    else
        path3 = [img_path,...
            num2str(curframe,'%04d'), '.bmp'];
    end
       
img3 = imread(path3); 
if ndims( img3 ) == 3
    img3 = rgb2gray( img3 );
end
curImg = double(img3);

tarImg = res_patch_ten(curframetartensor, curImg, patchSize, slideStep);
res_path = [des_path,num2str(curframindex, '%04d'), '.png'];
imwrite(tarImg / (max(tarImg(:)) + eps), res_path);
        if curframindex==initframe 
            for index=1:3
                tarImg = res_patch_ten(L(:,:,:,index), curImg, patchSize, slideStep);
                res_path = [des_path, num2str(index, '%04d'), '.png'];
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

toc
time=time+toc;
 %figure,imshow(tarImg,[])
end
res_base_path='D:\code\tensor_code\shiyan\benchmark_release(sea))_2\time_test_results\';
txt_path = [res_base_path, 'time_STT_easy',  num2str(tt) , '.txt'];
fid = fopen(txt_path, 'a');
times = 1.0 * time/44 ;
disp( num2str(times, '%.4f'));
fprintf(fid, '%s\n', num2str(times, '%.4f'),'times.');
fclose(fid);
end   

