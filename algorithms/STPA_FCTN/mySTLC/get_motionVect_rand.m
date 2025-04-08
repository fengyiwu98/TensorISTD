% *************************************************************************
% 这是一个采用穷尽搜索进行块匹配以计算背景运动矢量的函数
% 核心思想是随机选取 mbNum 个大小为 mbSize*mbSize 的宏块进行块匹配，以获取imgI相对于imgP的运动矢量
%（由imgP中核心块的位置指向在imgI中搜索到的与核心块相匹配的块的位置）
% 
% 输入：
% imgP - 当前帧
% imgI - 被搜索的帧
% mbSize - 宏块的边长
% mbNum - 宏块的数量
% p - 以当前宏块上下左右各p个像素为搜索区域，即搜索区域大小为：(2p+1)*(2p+1)
%
% 输出：
% motionVect - imgI相对于imgP的运动矢量（先行后列，行向下为正，列向右为正）
% *************************************************************************

function motionVect = get_motionVect_rand(imgP, imgI, mbSize, mbNum, p)

imgP = double(imgP);
imgI = double(imgI);

[row, col] = size(imgI);

valid_start_row = p + 1;
valid_end_row = row - p - mbSize;
valid_start_col = p + 1;
valid_end_col = col - p - mbSize;
% 随机 mbNum 个宏块的左上角坐标
rowPosArr = randi( [valid_start_row, valid_end_row], [mbNum, 1] );
colPosArr = randi( [valid_start_col, valid_end_col], [mbNum, 1] );

vectors = zeros(2, mbNum);  % 记录各宏块的运动矢量
costs = ones(2*p+1, 2*p+1) * 65537;  % 256 * 256 = 65536
% 2p+1行，2p+1列，初始值都为65537（一个大值）

mbCount = 1;  % 当前块在当前帧中的顺序编号
% (i, j)是当前宏块左上角的坐标
for idx = 1: mbNum
    
	i = rowPosArr(idx);
	j = colPosArr(idx);
    for m = -p : p        
        for n = -p : p
            refBlkVer = i + m;   % 参考块的垂直坐标
            refBlkHor = j + n;   % 参考块的水平坐标
            % 截取从当前帧和参考帧中分别截取当前块和参考块，计算MAD代价
            % MAD：每个像素点误差的绝对值之和 / 总像素个数
            costs(m+p+1,n+p+1) = costFuncMAD(imgP(i:i+mbSize-1, j:j+mbSize-1), ...
                imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
        end
    end
        
	% (i, j) - 当前块在当前帧中的绝对坐标
	% (p, p) - 以搜索区域左上角为原点，当前块在搜索区域中的相对坐标
	% (dx, dy) - 参考块在搜索区域中的相对坐标
	[dx, dy, ~] = minCost(costs); % 找出imgI中哪个宏块给了我们最小的代价
	% vector：当前块指向参考块的运动矢量
	vectors(1,mbCount) = dy-p-1;    % 向量的行坐标
	vectors(2,mbCount) = dx-p-1;    % 向量的列坐标
	mbCount = mbCount + 1;  % 当前块编号
	costs = ones(2*p + 1, 2*p +1) * 65537;  % 重新初始化代价
    
end

motionVect = round( mean(vectors,2) );

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Computes the Mean Absolute Difference (MAD) for the given two blocks
% Input
%       currentBlk : The block for which we are finding the MAD
%       refBlk : the block w.r.t. which the MAD is being computed
%       n : the side of the two square blocks
%
% Output
%       cost : The MAD for the two blocks
%
% Written by Aroh Barjatya
%
function cost = costFuncMAD(currentBlk,refBlk, n)
err = 0;
for i = 1:n
    for j = 1:n
        err = err + abs((currentBlk(i,j) - refBlk(i,j)));
    end
end
cost = err / (n*n);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Finds the indices of the cell that holds the minimum cost
%
% Input
%   costs : The matrix that contains the estimation costs for a macroblock
%
% Output
%   dx : the motion vector component in columns
%   dy : the motion vector component in rows
%
% Written by Aroh Barjatya

function [dx, dy, min] = minCost(costs)

[row, col] = size(costs);

% we check whether the current
% value of costs is less then the already present value in min. If its
% inded smaller then we swap the min value with the current one and note
% the indices.

min = 65537;

for i = 1:row
    for j = 1:col
        if (costs(i,j) < min)
            min = costs(i,j);
            dx = j; dy = i;
        end
    end
end

end
% -------------------------------------------------------------------------