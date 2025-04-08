function hog = setHoGParams(bins, rank, cellsz)
%SETHOGPARAMS 此处显示有关此函数的摘要
    % orign HoG struct is learned from standard by 'chrisjmccormick'.
    % bins - The number of bins to use in the histograms.
    % rank, nHC, nVC - 
    %   The number of cells horizontally and vertically.
    %   We force the sliding kernel as square with length of (rank*2+1)
    %   here.
    % cellsz - 
    %   Cell size in pixels (the cells are square), n*n pixels in each cell.
    % 
    % According to dual-layer sliding model given in section 3.2 in
    % Qian2020, we did adaptions based on original one:
    % T - 
    %   central cell of sliding block in which potential targets are
    %   scanned.
    % B - 
    %   non-central cell.
    % P - 
    %   diagonal directed shade blocks contains T alone its edge with diagonal
    %   length equals to ceil of halved diagonal length of whole block.
    %   Single contrast shade block composed of several cells.
    %   Contrast block here are compiled into 'INDEX-INDICATOR' within
    %   manner of cell index in local sampling window.
    
%   此处显示详细说明
%%
% Define the default HOG detector parameters (the size of the detection 
% window and the number of histogram bins).

% The number of bins to use in the histograms.
hog.numBins = bins;

% The number of cells horizontally and vertically.
% nHC * nVC cells in each single block.
nHC = (rank*2+1);
nVC = nHC;  % square block has same length and width
rkP = rank + 1;   % half diagonal length plus 1 equals to rank of contrast matrix
hog.numHorizCells = nHC;
hog.numVertCells = nVC;
% total member cell quantity in current HoG window
hog.ncell = hog.numHorizCells * hog.numVertCells; 

% Cell size in pixels (the cells are square).
hog.cellSize = cellsz;

% Compute the expected window size 
% (We removed 1 pixel border on all sides in original file).
hog.winR = (hog.numVertCells * hog.cellSize);
hog.winC = (hog.numHorizCells * hog.cellSize);
hog.winSize = [hog.winR, hog.winC];

% setup local sampling offset matrix
hog.Ty = floor(hog.winR/2)+1;
hog.Tx = floor(hog.winC/2)+1;
% sample offest matrix for y and x coords, y for Row and x for Col.
hog.smpofs = zeros(hog.winR, hog.winC, 2);  
for i = 1:hog.winR
    hog.smpofs(i,:,1) = i - hog.Ty;
end
for i = 1:hog.winC
    hog.smpofs(:,i,2) = i - hog.Tx;
end

% setup sampling section params
hog.nP = rkP^2; % Quantity of HoG contrast block equals to its total shade.
hog.block_cell_list = zeros(hog.nP);    % cell index list of contrast blocks

% allocate sampling cell index into contrast block
% from right-bottom cornor, we rotate block matrix alone clockwise order
% with target cell as center.

% setup cell index map
hog.cmp = zeros(nVC, nHC);
hog.cidx_pos = zeros(hog.ncell, 2); % position for single cell under block view 
cell_idx = 1;
for i = 1:nVC
    for j = 1:nHC
        hog.cmp(i,j) = cell_idx;
        hog.cidx_pos(cell_idx, 1) = i;
        hog.cidx_pos(cell_idx, 2) = j;
        cell_idx = cell_idx + 1;
    end
end

% cell position and cell index of central Target cell
hog.Ty_c = floor(nVC/2)+1;
hog.Tx_c = floor(nHC/2)+1;
hog.Tcidx = (hog.Ty_c - 1)*nHC + hog.Tx_c;

% regist sample cell index into contrast block list, with this registration
% list we can visit local histogram by cell via index indicator in contrast
% block.
blk_idx = 1;
for i = 0:rkP - 1
    for j = 0:rkP - 1
        % set left-top position of current contrast block 
        ofs_y = -hog.Ty_c + i;
        ofs_x = -hog.Tx_c + j;
        cell_idx = 1;
        for dr = 1:rkP
            for dc = 1:rkP
                % register a new member cell index from cell map into contrast 
                % list of current block 
                hog.block_cell_list(blk_idx,cell_idx) = hog.cmp(i + dr , j + dc);
                cell_idx = cell_idx + 1;
            end
        end
        blk_idx = blk_idx + 1;
    end
end

% init cell index map for grids in local sample matrix
hog.cellmap = zeros(hog.winR, hog.winC);
for i = 1:nVC
    for j = 1:nHC
        cell_idx = hog.cmp(i,j);    % pick cell index
        
        % fill cell-index and cell position into grids belong to current cell on cellmap
        for cell_y = 1:hog.cellSize
            l = (i-1)*hog.cellSize + cell_y;
            for cell_x = 1:hog.cellSize
                m = (j-1)*hog.cellSize + cell_x;
                hog.cellmap( l, m, 1) = cell_idx;
                hog.cellmap( l, m, 2) = cell_y;
                hog.cellmap( l, m, 3) = cell_x;
            end
        end
    end
end

end

