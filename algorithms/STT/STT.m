% This code is an implement of LogTFNN.
% Ref: H. -K. Liu, L. Zhang and H. Huang, "Small Target Detection inInfrared Videos Based on Spatio-Temporal Tensor Model", TGRS, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: 
% Email: 
% Copyright: University of Electronic Science and Technology of China
% Date: 2025/5/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef STT

    properties
        % para
        patchSize = 50;      
        slideStep = 30;       
         
        result;              % tar
    end
    

    
    methods
        %function obj = process(obj, inImg)
         %   img = inImg* 255;
         %   frame = 3;
         %   patchNum = 3; % patchNum X patchNum涓猵atch
         %   tem = fix(patchNum/2)
         %   len = length(img);
         %   img_num = len;
         %   initframe=4;
         %   endframe=img_num-3;
         
        function obj = process(obj, inseq_root_path, inidx3)
        frame = 3;       % 前后各frame帿
        patchNum = 3; 
        tem = fix(patchNum/2)
        %len = length(inImg);
        len = 100;
        img_num = len;
        initframe=4;
        endframe=img_num-3;
        
         path1 = [inseq_root_path,num2str(1,'%04d'), '.bmp']; 
%path1 = [inseq_root_path];
img1 = imread(path1); 

if ndims( img1 ) == 3
    img1 = rgb2gray( img1 );
end
            [imgHei, imgWid] = size(img1);

            rowPatchNum = ceil((imgHei - obj.patchSize) / obj.slideStep) + 1;
            colPatchNum = ceil((imgWid - obj.patchSize) / obj.slideStep) + 1;

            rowPosArr = [1 : obj.slideStep : (rowPatchNum - 1) * obj.slideStep, imgHei - obj.patchSize + 1];
            colPosArr = [1 : obj.slideStep : (colPatchNum - 1) * obj.slideStep, imgWid - obj.patchSize + 1];

            curcel = cell(rowPatchNum,colPatchNum,img_num);
            extend_curcel = cell(rowPatchNum+2*tem,colPatchNum+2*tem,img_num);

       
    %% 鎶妏atch鏀惧叆鍒癱ell涓?  
  
   for i = 1:img_num

              path = [inseq_root_path,num2str(i,'%04d'),  '.bmp'];

   % path = [file_path,num2str(i),'.png'];
    img2 = imread(path); 
   
    %%格式问题
    if ndims( img2 ) == 3
        img2 = rgb2gray( img2 );
    end
    img2 = double(img2);

    for row = 1:rowPatchNum
        for col = 1:colPatchNum
            tmp_patch = img2(rowPosArr(1,row) : rowPosArr(1,row) + obj.patchSize - 1, colPosArr(1,col) : colPosArr(1,col) + obj.patchSize - 1);
            curcel{row,col,i} = tmp_patch;
            extend_curcel{row+tem,col+tem,i} = tmp_patch;
        end
    end
  
    %鎵╁睍
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
  
%for curframindex = initframe:endframe 
 
 curframe =  inidx3+3;
%curframe =  curframindex;
tensorNum = patchNum*patchNum*(frame*2+1); 
curpatch = patchNum*patchNum*frame+patchNum*tem+tem+1;
curframetartensor = zeros(obj.patchSize, obj.patchSize, rowPatchNum*colPatchNum);
curframebartensor = zeros(obj.patchSize, obj.patchSize, rowPatchNum*colPatchNum);
t = 0;

for row = tem+1:rowPatchNum+tem  
    for col = tem+1:colPatchNum+tem
        
        k = 0;
        tmp_curpatchtensor = zeros(obj.patchSize, obj.patchSize, tensorNum);        
        %鏋勫缓浜嗗綋鍓嶇殑甯х殑褰撳墠鍧楃殑patch tensor
        for frameNum = curframe-frame:curframe+frame-3
            for patchRow = -tem:tem
                for patchCol = -tem:tem
                    k = k+1;
                    tmp_curpatchtensor(:,:,k) = extend_curcel{row+patchRow,col+patchCol,frameNum};
                end
            end           
        end
        
        %璋冩崲椤哄簭
        tmp_curpatchtensor1 = tmp_curpatchtensor(:,:,patchNum*patchNum*frame+1:patchNum*patchNum*(frame+1));
        tmp_curpatchtensor2 = tmp_curpatchtensor(:,:,patchNum*patchNum*(frame+1)+1:tensorNum);
        tmp_curpatchtensor(:,:,patchNum*patchNum*frame+1:patchNum*patchNum*2*frame) = tmp_curpatchtensor2;
        tmp_curpatchtensor(:,:,patchNum*patchNum*2*frame+1:tensorNum) = tmp_curpatchtensor1;

        lambda = 1/sqrt(patchNum*patchNum*(1+2*frame*obj.patchSize)+obj.slideStep);
        [tenBack,tenTar] = obj.TRPCA(tmp_curpatchtensor, lambda);  
       
        t = t+1;       
        curframetartensor(:,:,t) = tenTar(:,:,curpatch);
        curframebartensor(:,:,t) = tenBack(:,:,curpatch); 
        if inidx3==initframe 
            for index=1:3
                L(:,:,t,index)=  tenTar(:,:,(index-1)*patchNum*patchNum+patchNum*tem+tem+1);
            end
        end
         if inidx3==endframe 
            for index=1:3
                L(:,:,t,index)=  tenTar(:,:,(index+4-1)*patchNum*patchNum+patchNum*tem+tem+1);
            end
        end
             
    end
end

   path2 = [inseq_root_path,num2str(curframe,'%04d'),  '.bmp'];
curImg = imread(path2);   

%%格式问题
if ndims( curImg ) == 3
    curImg = rgb2gray( curImg );
end
curImg = double(curImg);

tarImg = obj.res_patch_ten(curframetartensor, curImg, obj.patchSize, obj.slideStep);
obj.result = tarImg;
        %end
        end

        function recImg = res_patch_ten(obj, patchTen, img, patchSize, slideStep)

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
        end

        function [L,S] = TRPCA(obj, X, lambda, opts)

    tol = 1e-1; %灏忎簬杩欎釜鏁板氨鍋滄
    max_iter = 20;
    rho = 1.05; 
    mu = 2*1e-3;
    max_mu = 1e10;
    DEBUG = 0;
    N = obj.rankN(X,0.1); %X鐨勭З


    if ~exist('opts', 'var')
        opts = [];
    end
    if isfield(opts, 'tol');         tol = opts.tol;              end
    if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
    if isfield(opts, 'rho');         rho = opts.rho;              end
    if isfield(opts, 'mu');          mu = opts.mu;                end
    if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
    if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
    if isfield(opts, 'N');           N = opts.N;                  end

    dim = size(X);
    L = zeros(dim);
    S = zeros(dim);
    Y = zeros(dim);

    for iter = 1 : max_iter

        preT = sum(S(:) > 0);
        
        % update L 鑳屾櫙
        R = -S+X-Y/mu;
        L = obj.prox_tnn(R,1/mu);

        % update S 鐩爣
        T = -L+X-Y/mu;
        S = obj.prox_l1(T, lambda/mu);
        
        % updata Y mu 鍏朵粬鍙傛暟
        dY = L+S-X;
        err = norm(dY(:))/norm(X(:));

        if DEBUG
            if iter == 1 || mod(iter, 1) == 0
                disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', err=' num2str(err)...
                    ',|T|0 = ' num2str(sum(S(:) > 0))]);
            end
        end

        currT = sum(S(:) > 0);

        if err < tol || (preT>0 && currT>0 && preT == currT)
            break;
        end
        
        Y = Y + dY*mu;
        mu = min(rho*mu,max_mu);
        
    end
end

function N = rankN(obj,X, ratioN)
    [~,~,n3] = size(X);
    D = obj.Unfold(X,n3,1); % 妯? 
    [~, S, ~] = svd(D, 'econ');
    [desS, ~] = sort(diag(S), 'descend');
    ratioVec = desS / desS(1);
    idxArr = find(ratioVec < ratioN);
    if idxArr(1) > 1
        N = idxArr(1) - 1;
    else
        N = 1;
    end
end
    
function [X] = Unfold(obj, X, dim, i )
     X = reshape(shiftdim(X,i-1), dim(i), []);
end

function x = prox_l1(obj, b,lambda)

x = max(0,b-lambda)+min(0,b+lambda);
x = max(x,0);

end



           function [X,tnn,trank] = prox_tnn(obj, Y,rho)


[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);
tnn = 0;
trank = 0;
        
% first frontal slice
[U,S,V] = svd(Y(:,:,1),'econ');
S = diag(S);
r = length(find(S>rho));
if r>=1
    S = S(1:r)-rho;
    X(:,:,1) = U(:,1:r)*diag(S)*V(:,1:r)';
    tnn = tnn+sum(S);
    trank = max(trank,r);
end
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    r = length(find(S>rho));
    if r>=1
        S = S(1:r)-rho;
        X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
        tnn = tnn+sum(S)*2;
        trank = max(trank,r);
    end
    X(:,:,n3+2-i) = conj(X(:,:,i));
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    r = length(find(S>rho));
    if r>=1
        S = S(1:r)-rho;
        X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
        tnn = tnn+sum(S);
        trank = max(trank,r);
    end
end
tnn = tnn/n3;
X = ifft(X,[],3);
           end
    end
end 