
function tenD = gen_patch(img, patchSize, slideStep)

if ~exist('patchSize', 'var')
    patchSize = 50;
end

if ~exist('slideStep', 'var')
    slideStep = 10;
end

    [imgHei, imgWid,n3] = size(img);
    
    rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
    colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;
    rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
    colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];
    patchTen = zeros(patchSize, patchSize, rowPatchNum * colPatchNum,n3);
for i=1:n3


    %% arrayfun version, identical to the following for-loop version
    %[meshCols, meshRows] = meshgrid(colPosArr, rowPosArr);
    %idx_fun = @(row,col) img(row : row + patchSize - 1, col : col + patchSize - 1);
    %patchCell = arrayfun(idx_fun, meshRows, meshCols, 'UniformOutput', false);
    %patchTen = cat(3, patchCell{:});

    %% for-loop version
     
     k = 0;
     for row = rowPosArr
         for col = colPosArr
             k = k + 1;
             tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1,i);
             patchTen(:, :, k,i) = tmp_patch;
         end
     end
end
[~,~,n2,~] = size(patchTen);
tenD=zeros(patchSize,patchSize,27,(colPatchNum-2)*(rowPatchNum-2));
ii=1;
square_arr=[];
for i=1:rowPatchNum-2
    square_arr1=i*colPatchNum+2:(i+1)*colPatchNum-1;
    square_arr=[square_arr square_arr1 ];
end
for j=square_arr
    %第一张
    tenD3(:,:,1)=patchTen(:,:,j-colPatchNum-1,1);
    tenD3(:,:,2)=patchTen(:,:,j-colPatchNum,1);
    tenD3(:,:,3)=patchTen(:,:,j-colPatchNum+1,1);
    tenD3(:,:,4)=patchTen(:,:,j-1,1);
    tenD3(:,:,5)=patchTen(:,:,j,1);
    tenD3(:,:,6)=patchTen(:,:,j+1,1);
    tenD3(:,:,7)=patchTen(:,:,j+colPatchNum-1,1);
    tenD3(:,:,8)=patchTen(:,:,j+colPatchNum,1);
    tenD3(:,:,9)=patchTen(:,:,j+colPatchNum+1,1);
    %第二张
    tenD3(:,:,10)=patchTen(:,:,j-colPatchNum-1,2);
    tenD3(:,:,11)=patchTen(:,:,j-colPatchNum,2);
    tenD3(:,:,12)=patchTen(:,:,j-colPatchNum+1,2);
    tenD3(:,:,13)=patchTen(:,:,j-1,2);
    tenD3(:,:,14)=patchTen(:,:,j,2);
    tenD3(:,:,15)=patchTen(:,:,j+1,2);
    tenD3(:,:,16)=patchTen(:,:,j+colPatchNum-1,2);
    tenD3(:,:,17)=patchTen(:,:,j+colPatchNum,2);
    tenD3(:,:,18)=patchTen(:,:,j+colPatchNum+1,2);
    %第三张
    tenD3(:,:,19)=patchTen(:,:,j-colPatchNum-1,3);
    tenD3(:,:,20)=patchTen(:,:,j-colPatchNum,3);
    tenD3(:,:,21)=patchTen(:,:,j-colPatchNum+1,3);
    tenD3(:,:,22)=patchTen(:,:,j-1,3);
    tenD3(:,:,23)=patchTen(:,:,j,3);
    tenD3(:,:,24)=patchTen(:,:,j+1,3);
    tenD3(:,:,25)=patchTen(:,:,j+colPatchNum-1,3);
    tenD3(:,:,26)=patchTen(:,:,j+colPatchNum,3);
    tenD3(:,:,27)=patchTen(:,:,j+colPatchNum+1,3);
    %第四张
    tenD3(:,:,28)=patchTen(:,:,j-colPatchNum-1,4);
    tenD3(:,:,29)=patchTen(:,:,j-colPatchNum,4);
    tenD3(:,:,30)=patchTen(:,:,j-colPatchNum+1,4);
    tenD3(:,:,31)=patchTen(:,:,j-1,4);
    tenD3(:,:,32)=patchTen(:,:,j,4);
    tenD3(:,:,33)=patchTen(:,:,j+1,4);
    tenD3(:,:,34)=patchTen(:,:,j+colPatchNum-1,4);
    tenD3(:,:,35)=patchTen(:,:,j+colPatchNum,4);
    tenD3(:,:,36)=patchTen(:,:,j+colPatchNum+1,4);
    %第五张
    tenD3(:,:,37)=patchTen(:,:,j-colPatchNum-1,5);
    tenD3(:,:,38)=patchTen(:,:,j-colPatchNum,5);
    tenD3(:,:,39)=patchTen(:,:,j-colPatchNum+1,5);
    tenD3(:,:,40)=patchTen(:,:,j-1,5);
    tenD3(:,:,41)=patchTen(:,:,j,5);
    tenD3(:,:,42)=patchTen(:,:,j+1,5);
    tenD3(:,:,43)=patchTen(:,:,j+colPatchNum-1,5);
    tenD3(:,:,44)=patchTen(:,:,j+colPatchNum,5);
    tenD3(:,:,45)=patchTen(:,:,j+colPatchNum+1,5);
    tenD(:,:,1:45,ii)=tenD3;
    
    ii=ii+1;
end

     


