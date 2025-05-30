 % This code is an implement of STPA_FCTN.
% Ref: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: 
% Email: 
% Copyright:  University of Electronic Science and Technology of China
% Date: 2025/5/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef STPA_FCTN

    properties
        patchSize =70;  
        slideStep =70;
        timeSize =15;  

        result;              % tar
    end

    methods
        function obj = process(obj, inseq_root_path, inpng_res_base_path, inlen)
            tol = [1e-3];
            ext_name='.bmp';
            start_idx = 1;
            %len = 100;
            i_num = ceil(inlen/obj.timeSize);
            for i=1:i_num
                for j=1:obj.timeSize
                    if i>floor(inlen/obj.timeSize)
                        cur_idx = j + inlen - obj.timeSize + start_idx - 1;
                    else
                        cur_idx = (i - 1) * obj.timeSize + j + start_idx - 1;
                    end

                    img = imread( [inseq_root_path num2str(cur_idx,'%04d') '.bmp'] );
                    if ndims(img) == 3
                        img = rgb2gray(img);
                    end
                    [m,n]=size(img);

                    %prior image
                    %STPA algorithm
                    last_i = i_num;
                    last_j = obj.timeSize - 2;
                    if (i==1 && j<3)||(i==last_i && j>last_j)
                        prior_map=obj.mySLC(inseq_root_path,ext_name,cur_idx);
                    else
                        prior_map=obj.mySTLC(inseq_root_path,ext_name,cur_idx);
                    end

                    %if (i==1 && j<3)||(i==last_i && j>last_j)
                    %prior_map=mySLC(img_path,ext_name,cur_idx);
                    %else
                    %prior_map=mySTLC(img_path,ext_name,cur_idx);
                    %end
                    %直接读取结果图像
                    %prior_map=imread([prior_result_path num2str(cur_idx,'%04d') '.png']);

                    % 4D tensor
                    prior_map=im2double(prior_map);
                    img = im2double(img);
                    D = obj.gen_patch_ten(img, obj.patchSize, obj.slideStep);
                    prior_3D=obj.gen_patch_ten(prior_map, obj.patchSize, obj.slideStep);
                    img_4D(:,:,j,:)=D;
                    prior_4D(:,:,j,:)=prior_3D;
                end

                
                [tenBack,tenTar] = obj.FCTN(img_4D,prior_4D, tol);
                %[tenBack,tenTar] = FCTN(img_4D,prior_4D); %prior weight
                %[tenBack,tenTar] = FCTN(img_4D,tol);%without prior weight
              

                for j = 1: obj.timeSize
                    if i>floor(inlen/obj.timeSize)
                        cur_idx = inlen +j- obj.timeSize;
                    else
                        cur_idx = (i-1) * obj.timeSize + j + start_idx - 1;
                    end
                    tarImg = obj.res_patch_ten_mean (tenTar(:,:,j,:), img, obj.patchSize, obj.slideStep);
                    T = tarImg / (max(tarImg(:)) + eps);

                    res_path = [inpng_res_base_path, num2str(cur_idx,'%04d'), '.png'];
                    imwrite(T, res_path);
                end

            end
        end

        function Smap = mySTLC(obj, file_path, ext_name, img_idx)
            % tic
            SLC = obj.mySLC(file_path, ext_name, img_idx);
            TLC = obj.myTLC2(file_path, ext_name, img_idx);
            Smap = SLC .* TLC;
            Smap = ( Smap - min(Smap(:)) ) / (  max(Smap(:)) -  min(Smap(:)) );
            %Smap = Smap/ (  max(Smap(:)) -  min(Smap(:)) );
            imshow(Smap);
        end

        function SLC = mySLC(obj, file_path, ext_name, img_idx)

            path = fullfile([file_path, num2str(img_idx,'%04d'), ext_name]);
            img = imread(path);
            if ndims( img ) == 3
                img = rgb2gray(img);
            end
            img = double(img);

            x = [-1, 0, 1];

            sigma_o = 1;
            y_o = (1./(sqrt(2.*pi).*sigma_o)) .* exp(-x.^2./(2.*sigma_o.^2));
            w_o = y_o ./ sum(y_o);

            sigma_i = 0.3;
            y_i = (1./(sqrt(2.*pi).*sigma_i)) .* exp(-x.^2./(2.*sigma_i.^2));
            w_i = y_i ./ sum(y_i);

            m1 = zeros(9, 9);
            m1(1,1)=-w_o(1); m1(2,2)=-w_o(2); m1(3,3)=-w_o(3);
            m1(4,4)=w_i(1); m1(5,5)=w_i(2); m1(6,6)=w_i(3);

            m2 = zeros(9, 9);
            m2(1,5)=-w_o(1); m2(2,5)=-w_o(2); m2(3,5)=-w_o(3);
            m2(4,5)=w_i(1); m2(5,5)=w_i(2); m2(6,5)=w_i(3);

            m3 = zeros(9, 9);
            m3(1,9)=-w_o(1); m3(2,8)=-w_o(2); m3(3,7)=-w_o(3);
            m3(4,6)=w_i(1); m3(5,5)=w_i(2); m3(7,4)=w_i(3);

            m4 = zeros(9, 9);
            m4(5,9)=-w_o(1); m4(5,8)=-w_o(2); m4(5,7)=-w_o(3);
            m4(5,6)=w_i(1); m4(5,5)=w_i(2); m4(5,4)=w_i(3);

            m5 = zeros(9, 9);
            m5(9,9)=-w_o(1); m5(8,8)=-w_o(2); m5(7,7)=-w_o(3);
            m5(6,6)=w_i(1); m5(5,5)=w_i(2); m5(4,4)=w_i(3);

            m6 = zeros(9, 9);
            m6(9,5)=-w_o(1); m6(8,5)=-w_o(2); m6(7,5)=-w_o(3);
            m6(6,5)=w_i(1); m6(5,5)=w_i(2); m6(4,5)=w_i(3);

            m7 = zeros(9, 9);
            m7(9,1)=-w_o(1); m7(8,2)=-w_o(2); m7(7,3)=-w_o(3);
            m7(6,4)=w_i(1); m7(5,5)=w_i(2); m7(4,6)=w_i(3);

            m8 = zeros(9, 9);
            m8(5,1)=-w_o(1); m8(5,2)=-w_o(2); m8(5,3)=-w_o(3);
            m8(5,4)=w_i(1); m8(5,5)=w_i(2); m8(5,6)=w_i(3);

            S1 = abs(imfilter(img, m1, 'replicate'));
            S2 = abs(imfilter(img, m2, 'replicate'));
            S3 = abs(imfilter(img, m3, 'replicate'));
            S4 = abs(imfilter(img, m4, 'replicate'));
            S5 = abs(imfilter(img, m5, 'replicate'));
            S6 = abs(imfilter(img, m6, 'replicate'));
            S7 = abs(imfilter(img, m7, 'replicate'));
            S8 = abs(imfilter(img, m8, 'replicate'));

            S = min( min( min(S1,S2), min(S3,S4) ), min( min(S5,S6), min(S7,S8) ) );
            SLC = S / max(S(:));


        end


        function TLC = myTLC2(obj, file_path, ext_name, img_idx)

            path_cur = fullfile([file_path, num2str(img_idx,'%04d'), ext_name]);
            img_cur = imread(path_cur);
            if ndims( img_cur ) == 3
                img_cur = double( rgb2gray(img_cur) );
            end
            [imgR, imgC] = size(img_cur);

            frame_data = zeros(imgR, imgC, 4);
            motionVect = zeros(2, 4);

            for k = -2: 2

                temp = imread(fullfile([file_path, num2str(img_idx+k,'%04d'), ext_name]));
                if ndims(temp) == 3
                    temp = rgb2gray(temp);
                end

                if k < 0
                    frame_data(:,:,k+3) = double(temp);
                    motionVect(:,k+3) = obj.get_motionVect_rand2(img_cur, temp, 64, 3, -8*k);
                end
                if k > 0
                    frame_data(:,:,k+2) = double(temp);
                    motionVect(:,k+2) = obj.get_motionVect_rand2(img_cur, temp, 64, 3, 8*k);
                end

            end

            S = zeros(imgR, imgC);
            %%%%% 根据运动向量获取坐标？
            for y = 1:imgR
            	for x = 1:imgC
                    y1 = min( max(y+motionVect(1,1), 1), imgR);
                    y2 = min( max(y+motionVect(1,2), 1), imgR);
                    y3 = min( max(y+motionVect(1,3), 1), imgR);
                    y4 = min( max(y+motionVect(1,4), 1), imgR);
                    x1 = min( max(x+motionVect(2,1), 1), imgC);
                    x2 = min( max(x+motionVect(2,2), 1), imgC);
                    x3 = min( max(x+motionVect(2,3), 1), imgC);
                    x4 = min( max(x+motionVect(2,4), 1), imgC);
                    mt = 0.25 * ( frame_data(y1,x1,1) + frame_data(y4,x4,4) ) / 2  + ...
                        0.75 * ( frame_data(y2,x2,2) + frame_data(y3,x3,3) ) / 2;
                    % S(y,x) = abs( img_cur(y,x) - mt );
                    S(y,x) = max( img_cur(y,x) - mt,0 );
            	end
            end

            TLC = S  / max(S(:));

            %imshow(TLC);

        end

        function motionVect = get_motionVect_rand2(obj, imgP, imgI, mbSize, mbNum, p)

            imgP = double(imgP);
            imgI = double(imgI);

            [row, col] = size(imgI);

            valid_start_row = p + 1;
            valid_end_row = row - p - mbSize;
            valid_start_col = p + 1;
            valid_end_col = col - p - mbSize;
            % 固定 mbNum 个宏块的左上角坐标
            rowNum = floor(sqrt(mbNum));%横向选取sqrt()个
            colNum = ceil(mbNum/rowNum);
            rowSlide = fix((valid_end_row - valid_start_row)/rowNum); %横向宏块点间隔
            colSlide = fix((valid_end_col - valid_start_col)/colNum);%纵向间隔
            rowPosArr = repmat([valid_start_row : rowSlide : rowSlide * rowNum],[1,colNum]);
            colPosArr = repelem([valid_start_col : colSlide : colSlide * colNum],rowNum);
            %rowPosArr = randi( [valid_start_row, valid_end_row], [mbNum, 1] );
            %colPosArr = randi( [valid_start_col, valid_end_col], [mbNum, 1] );

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
                        costs(m+p+1,n+p+1) = obj.costFuncMAD(imgP(i:i+mbSize-1, j:j+mbSize-1), ...
                            imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                    end
                end

            	% (i, j) - 当前块在当前帧中的绝对坐标
            	% (p, p) - 以搜索区域左上角为原点，当前块在搜索区域中的相对坐标
            	% (dx, dy) - 参考块在搜索区域中的相对坐标
            	[dx, dy, ~] = obj.minCost(costs); % 找出imgI中哪个宏块给了我们最小的代价
            	% vector：当前块指向参考块的运动矢量
            	vectors(1,mbCount) = dy-p-1;    % 向量的行坐标
            	vectors(2,mbCount) = dx-p-1;    % 向量的列坐标
            	mbCount = mbCount + 1;  % 当前块编号
            	costs = ones(2*p + 1, 2*p +1) * 65537;  % 重新初始化代价
                %for i = 1:mbNum
                %   fprintf('paches %d: x=%d, y=%d\n', i, rowPosArr(i), colPosArr(i));
                %end
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
        function cost = costFuncMAD(obj,currentBlk,refBlk, n)
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

        function [dx, dy, min] = minCost(obj,costs)

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

        function patchTen = gen_patch_ten(obj, img, patchSize, slideStep)

            % 2017-07-31
            % This matlab code implements the RIPT model for infrared target-background
            % separation.
            %
            % Yimian Dai. Questions? yimian.dai@gmail.com
            % Copyright: College of Electronic and Information Engineering,
            %            Nanjing University of Aeronautics and Astronautics


            if ~exist('patchSize', 'var')
                patchSize = 50;
            end

            if ~exist('slideStep', 'var')
                slideStep = 15;
            end

            % img = reshape(1:9, [3 3])
            % img = reshape(1:12, [3 4])
            % patchSize = 2;
            % slideStep = 1;
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

        end

        %function [B,T,RC]=FCTN(D,W)
        function [B,T,RC]=FCTN(obj, D, W, tol)
            %function [B,T,RC]=FCTN(D,tol)
            %% initialize parameters

            H =0.015; %0.005 0.01 0.015 0.02 0.025 0.03
            %tol=1e-2;  
            maxiter=200;
            lambda1=0;
            %lambda2 = 10000;

            epsilon=1e-3;
            miu = 1e-4; %gamma %beta和gamma都设置成了miu;

            SR=1;
            Nway=size(D);
            N=ndims(D);
            L=nchoosek(N,2)/2; %取组合数之后再除2,或L直接设置为3；
            sk=zeros(L,1);
            Wmap=ones(Nway);
            RC=nan(maxiter,2);
            weight=1./L*ones(1,L); %delta


            %% initialization

            B=rand(Nway);
            T=zeros(Nway);
            E=zeros(Nway);
            %No=zeros(Nway);
            U=cell(L,1);
            C=cell(L,1);
            orderL=[1 2 3 4;3 1 2 4;2 3 1 4];

            for n=1:L
                U{n}=B;
                C{n}=zeros(Nway);
                order=orderL(n,:);
                lambda1=lambda1+H/sqrt(SR*max([prod(Nway(order(1:2))),prod(Nway(order(3:N)))]));
                %sk(n)=min([prod(Nway(order(1:L))),prod(Nway(order(L+1:N)))]);%添加？？
            end


            %% ADMM algorithm
            t=cputime;
            for i=1:maxiter
                preT=sum(T(:)>0);

                % update u^(n) (auxiliary variables of low-rank part)
                u_cs=zeros(Nway);
                z_cs=zeros(Nway);

                for n=1:L
                    order=orderL(n,:);
                    m=permute(B-C{n}/miu,order);
                    M=reshape(m,prod(Nway(order(1:2))),[]);
                    [M,sk(n)]=obj.shrink_matrix(M,weight(n)/miu,sk(n),epsilon,false);
                    m=reshape(M,Nway(order));
                    U{n}=ipermute(m,order);
                    u_cs=u_cs+U{n};
                    z_cs=z_cs+C{n};
                end

                % update B (low-rank part)
                B=(u_cs+z_cs/miu+(D-T-E/miu))./(L+1);
                u_cs=zeros(Nway);
                z_cs=zeros(Nway);

                % update T (sparse part)
                T=obj.shrink_vector((D-B-E/miu),Wmap.*lambda1/miu);
                %T=shrink_vector((D-B-E/miu),lambda1/miu); %without prior weight

                %update W(weignt part)
                Wmap = min( 1./(abs(T)+0.01)./W, 100);


                % update C^(n)
                for n=1:L
                    C{n}=C{n}+miu*(U{n}-B);
                end

                % update E
                E=E+miu*(B+T-D);

                %update N
                %No=miu*(D-B-T+E/miu)/(miu+2*lambda2);

                % evaluate recovery accuracy
                miu = min(miu*1.2,1e10);
                currT = sum(T(:) > 0);
                errList = norm(D(:)-B(:)-T(:)) / norm(D(:));
                %fprintf('tensor ring iterations = %d   difference=%f\n', i, errList);
                disp(['iter ' num2str(i)  ...
                    ', err=' num2str(errList)...
                    ',|T|0 = ' num2str(currT)]);

                if errList < tol ||(preT>0 && currT>0 && preT == currT)
                    break;
                end

            end

        end

        function [M1,num1]=shrink_matrix(obj,M0,tau,num0,epsilon,flag)
            [m,n]=size(M0);
            if flag
                if m<1.5*n
                    [U,Sigma2]=eigs(M0*M0',num0,'largestreal','MaxIterations',500);
                    sigma=sqrt(diag(Sigma2));
                    s=sigma-tau;
                    idx_s=find(s>0,1,'last');
                    num1=find(s>=epsilon*s(1),1,'last');
                    S_shrink=diag(s(1:idx_s)./sigma(1:idx_s));
                    M1=U(:,1:idx_s)*S_shrink*U(:,1:idx_s)'*M0;
                    return
                end
                if m>1.5*n
                    [M1,num1]=shrink_matrix(M0',tau,num0,epsilon,flag);
                    M1=M1';
                    return
                end
                [U,S,V]=svds(M0,num0,'largest','MaxIterations',500);
                s=diag(S)-tau;
                idx_s=find(s>0,1,'last');
                num1=find(s>=epsilon*s(1),1,'last');
                S_shrink=diag(s(1:idx_s));
                M1=U(:,1:idx_s)*S_shrink*V(:,1:idx_s)';
            else
                [U,S,V]=svd(M0,'econ');
                s=diag(S);
                idx=find(s>tau,1,'last');
                M1=U(:,1:idx)*diag(s(1:idx)-tau)*V(:,1:idx)';
                num1=num0;
            end
        end

        function v1=shrink_vector(obj,v0,tau)
            v1=sign(v0).*max(abs(v0)-tau,0);
        end

        function recImg = res_patch_ten_mean(obj, patchTen, img, patchSize, slideStep)

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
            for col = colPosArr
                for row = rowPosArr
                    k = k + 1;
                    tmpPatch = reshape(patchTen(:, :, k), [patchSize, patchSize]);
                    accImg(row : row + patchSize - 1, col : col + patchSize - 1) = tmpPatch;
                    weiImg(row : row + patchSize - 1, col : col + patchSize - 1) = onesMat;
                end
            end

            recImg = accImg ./ weiImg;
        end
    end
end