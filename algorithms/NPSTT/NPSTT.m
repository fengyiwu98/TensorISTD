 % This code is an implement of NPSTT.
% Ref: G. Wang, B. Tao, X. Kong and Z. Peng., "Infrared Small TargetDetection Using Nonoverlapping Patch Spatial–Temporal Tensor Factorization With Capped Nuclear Norm Regularization.",TGRS, 2022.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: 
% Email: 
% Copyright:  University of Electronic Science and Technology of China
% Date: 2025/5/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef NPSTT

    properties
        % para
        patchSize = 70;
        slideStep = 70;

        result;              % tar
    end

    methods

        function obj = process(obj, inseq_root_path, inpng_res_base_path, inlen)
            frame = 3;       % 前后各frame
            patchNum = 3;
            tem = fix(patchNum/2);

            img_num = inlen;
            initframe=4;
            endframe=img_num-3;

            path1 = [inseq_root_path,num2str(1,'%04d'), '.bmp'];
            img1 = imread(path1);

            if ndims( img1 ) == 3
                img1 = rgb2gray( img1 );
            end
            [imgHei, imgWid] = size(img1);
            rowPatchNum = ceil((imgHei - obj.patchSize) / obj.slideStep) + 1;  % patch行数
            colPatchNum = ceil((imgWid - obj.patchSize) / obj.slideStep) + 1;  % patch列数

            rowPosArr = [1 : obj.slideStep : (rowPatchNum - 1) * obj.slideStep, imgHei - obj.patchSize + 1];
            colPosArr = [1 : obj.slideStep : (colPatchNum - 1) * obj.slideStep, imgWid - obj.patchSize + 1];
            % 创建存放图片的Patch的cell
            % extend_curcel 为了解决边界问题
            curcel = cell(rowPatchNum,colPatchNum,img_num);  % 行数*列数*帧数
            extend_curcel = cell(rowPatchNum+2*tem,colPatchNum+2*tem,img_num);  % padding

            for i = 1:img_num
                path = [inseq_root_path,num2str(i,'%04d'),  '.bmp'];
                img2 = imread(path);

                %%格式问题
                if ndims( img2 ) == 3
                    img2 = rgb2gray( img2 );
                end
                img2 = double(img2);


                for row = 1:rowPatchNum
                    for col = 1:colPatchNum
                        % 在第i张图像中截取第row行第col列的patch
                        tmp_patch = img2(rowPosArr(1,row) : rowPosArr(1,row) + obj.patchSize - 1, colPosArr(1,col) : colPosArr(1,col) + obj.patchSize - 1);
                        curcel{row,col,i} = tmp_patch;
                        extend_curcel{row+tem,col+tem,i} = tmp_patch;
                    end
                end

                % 扩展（复制相邻像素）
                % 四角
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
                % 四边
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

            %for curframindex=initframe:endframe
            for curframindex=initframe:endframe
                curframe =  curframindex;
                %curframe =  inidx3+3;
                %% 构建NPSTT
                tensorNum = patchNum*patchNum*(frame*2+1);  % 7帧
                % 当前块在63个块中的位置
                % curpatch = patchNum*patchNum*2*frame+patchNum*tem+tem+1;  %错了？？＿
                curpatch = patchNum*patchNum*frame+patchNum*tem+tem+1;
                % 目标与背景张量的大小均与当前帧所有patch堆叠成的张量大小丿??
                curframetartensor = zeros(obj.patchSize, obj.patchSize, rowPatchNum*colPatchNum);
                curframebartensor = zeros(obj.patchSize, obj.patchSize, rowPatchNum*colPatchNum);
                t = 0;

                for row = tem+1:rowPatchNum+tem  % 选定当前块（不迨??padding部分＿
                    for col = tem+1:colPatchNum+tem

                        % 获取对应当前的帧的当前块的tensor
                        k = 0;
                        tmp_curpatchtensor = zeros(obj.patchSize, obj.patchSize, tensorNum);
                        for frameNum = curframe-frame:curframe+frame-3  % 与当前帧相关的前后帧
                            for patchRow = -tem:tem  % 当前块及其空间周围块
                                for patchCol = -tem:tem
                                    k = k+1;
                                    tmp_curpatchtensor(:,:,k) = extend_curcel{row+patchRow,col+patchCol,frameNum};
                                end
                            end
                        end

                        lambda = 0.5/sqrt(patchNum*patchNum*(1+2*frame)*obj.patchSize);
                        [tenBack,tenTar] = obj.tCSVT(tmp_curpatchtensor,lambda);
                        t = t+1;
                        curframetartensor(:,:,t) = tenTar(:,:,curpatch);
                        curframebartensor(:,:,t) = tenBack(:,:,curpatch);
                        %if curframindex==initframe
                        if curframindex==initframe
                            for index=1:3
                                L(:,:,t,index)=  tenTar(:,:,(index-1)*patchNum*patchNum+patchNum*tem+tem+1);
                            end
                        end
                        %if curframindex==endframe
                        if curframindex==endframe
                            for index=1:3
                                L(:,:,t,index)=  tenTar(:,:,(index+4-1)*patchNum*patchNum+patchNum*tem+tem+1);
                            end
                        end
                    end
                end


                path2 = [inseq_root_path,...
                    num2str(curframe,'%04d'),  '.bmp'];
                %path2 = [file_path,num2str(curframe),'.png'];
                curImg = imread(path2);

                %%格式问题
                if ndims( curImg ) == 3
                    curImg = rgb2gray( curImg );
                end
                curImg = double(curImg);

                tarImg = obj.res_patch_ten(curframetartensor, curImg, obj.patchSize, obj.slideStep);
                %backImg = res_patch_ten(curframebartensor, curImg, patchSize, slideStep);
                res_path = [inpng_res_base_path, num2str(curframe, '%04d'), '.png'];
                imwrite(tarImg / (max(tarImg(:)) + eps), res_path);

                if curframindex==initframe
                    for index=1:3
                        tarImg = obj.res_patch_ten(L(:,:,:,index), curImg, obj.patchSize, obj.slideStep);
                        res_path = [inpng_res_base_path,num2str(index, '%04d'), '.png'];
                        imwrite(tarImg / (max(tarImg(:)) + eps), res_path);
                    end
                end
                if curframindex==endframe
                    for index=1:3
                        tarImg = obj.res_patch_ten(L(:,:,:,index), curImg, obj.patchSize, obj.slideStep);
                        res_path = [inpng_res_base_path, num2str(index+endframe, '%04d'), '.png'];
                        imwrite(tarImg / (max(tarImg(:)) + eps), res_path);
                    end
                end

            end
        end

        function [L,S] = tCSVT(obj, X,lambda)

            tol = 1e-2;
            max_iter = 500;  
            rho = 1.05;
            mu = 1e-3;
            threshold = 3;
            max_mu = 1e10;
            DEBUG = 0;

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

            N = ones(1,1+round(dim(1,3)/2));  

            for iter = 1 : max_iter

                preT = sum(S(:) > 0);  

                % update L 
                R =  -S+X-Y/mu;
                L = obj.t_CSVT(R,N,mu);  
                for index = 1:1+round(dim(1,3)/2)
                   
                    [~,S,~] = svd(L(:,:,index),'econ');
                    diagS = diag(S);
                    thresh = max(diagS-threshold, 0);
                    
                    N(1,index) = size(find(thresh>0),1);
                end

                % update S 
                T = -L+X-Y/mu;
                S = obj.prox_l1(T, lambda/mu);

                % updata Y
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



        function x = prox_l1(obj,b,lambda)

            x = max(0,b-lambda)+min(0,b+lambda);
            x = max(x,0);

        end

        function [X] = t_CSVT(obj,Y,N,mu)
            [n1,n2,n3] = size(Y);
            X = zeros(n1,n2,n3);
            Y = fft(Y,[],3);
            tau = 1/mu;


            [U,S,V] = svd(Y(:,:,1),'econ');
            diagS = diag(S);
            [desS, sIdx] = sort(diagS, 'descend');
            [desU, desV] = deal(U(:, sIdx), V(:, sIdx));
            [U1, diagS1, V1] = deal(desU(:, 1:N(1,1)), desS(1:N(1,1)), desV(:, 1:N(1,1)));
            [U2, diagS2, V2] = deal(desU(:, N(1,1)+1:end), desS(N(1,1)+1:end), desV(:, N(1,1)+1:end));
            threshS2 = max(diagS2-tau, 0);
            X(:,:,1) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';


            % i=2,...,halfn3
            halfn3 = round(n3/2);
            for i = 2 : halfn3
                [U,S,V] = svd(Y(:,:,i),'econ');
                diagS = diag(S);
                [desS, sIdx] = sort(diagS, 'descend');
                [desU, desV] = deal(U(:, sIdx), V(:, sIdx));
                [U1, diagS1, V1] = deal(desU(:, 1:N(1,i)), desS(1:N(1,i)), desV(:, 1:N(1,i)));
                [U2, diagS2, V2] = deal(desU(:, N(1,i)+1:end), desS(N(1,i)+1:end), desV(:, N(1,i)+1:end));
                threshS2 = max(diagS2-tau, 0);
                X(:,:,i) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';
                X(:,:,n3+2-i) = conj(X(:,:,i));
            end

            % if n3 is even
            if mod(n3,2) == 0
                i = halfn3+1;
                [U,S,V] = svd(Y(:,:,i),'econ');
                diagS = diag(S);
                [desS, sIdx] = sort(diagS, 'descend');
                [desU, desV] = deal(U(:, sIdx), V(:, sIdx));
                [U1, diagS1, V1] = deal(desU(:, 1:N(1,2)), desS(1:N(1,2)), desV(:, 1:N(1,2)));
                [U2, diagS2, V2] = deal(desU(:, N(1,2)+1:end), desS(N(1,2)+1:end), desV(:, N(1,2)+1:end));
                threshS2 = max(diagS2-tau, 0);
                X(:,:,i) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';
            end

            X = ifft(X,[],3);
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

    end
end
