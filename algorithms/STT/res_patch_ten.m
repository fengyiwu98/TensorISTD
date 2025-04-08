 function recImg = res_patch_ten(patchTen, img, patchSize, slideStep)

% IPT复原模型

[imgHei, imgWid] = size(img);

rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];

%% for-loop version
accImg = zeros(imgHei, imgWid);
weiImg = zeros(imgHei, imgWid);

k = 0;
onesMat = ones(patchSize, patchSize);
for row = rowPosArr
    for  col = colPosArr
        k = k + 1;
        tmpPatch = reshape(patchTen(:, :, k), [patchSize, patchSize]);
%         accImg(row : row + patchSize - 1, col : col + patchSize - 1) = tmpPatch;
%         weiImg(row : row + patchSize - 1, col : col + patchSize - 1) = onesMat;
        accImg(row : row + patchSize - 1, col : col + patchSize - 1) = tmpPatch + accImg(row : row + patchSize - 1, col : col + patchSize - 1);
        weiImg(row : row + patchSize - 1, col : col + patchSize - 1) = onesMat +  weiImg(row : row + patchSize - 1, col : col + patchSize - 1);
    end
end

recImg = accImg ./ weiImg;
