function patchTen = gen_patch_ten(img, patchSize, slideStep)

if ~exist('patchSize', 'var')
    patchSize = 50;
end

if ~exist('slideStep', 'var')
    slideStep = 15;
end

[imgHei, imgWid] = size(img);

rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];

%% arrayfun version, identical to the following for-loop version
[meshCols, meshRows] = meshgrid(colPosArr, rowPosArr);
idx_fun = @(row,col) img(row : row + patchSize - 1, col : col + patchSize - 1);
patchCell = arrayfun(idx_fun, meshRows, meshCols, 'UniformOutput', false);
patchTen = cat(3, patchCell{:});

