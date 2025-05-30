 % This code is an implement of TR.
% Ref: 
% F. Wu, H. Yu, A. Liu, J. Luo, and Z. Peng, “Infrared small target detection using spatiotemporal 4-d tensor train and ring unfolding,” IEEE Trans. Geosci. Remote Sens., vol. 61, pp. 1–22, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: 
% Email: 
% Copyright:  University of Electronic Science and Technology of China
% Date: 2025/5/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef TR

    properties
        patchSize = 70;
        slideStep = 70;
        start_idx = 1;
        %length = 100;
        L = 15;

        result;              % tar
    end

    methods
        function obj = process(obj, inseq_root_path, inpng_res_base_path, inlen)
            for i = 1:floor(inlen/obj.L)
                for j=1:obj.L
                    cur_idx = (i-1) * obj.L + j + obj.start_idx - 1;
                    img = imread( [inseq_root_path num2str(cur_idx,'%04d') '.bmp'] );
                    if ndims(img) == 3
                        img = rgb2gray(img);
                    end
                    img = im2double(img);
                    D = obj.gen_patch_ten(img, obj.patchSize, obj.slideStep);
                    tensor_4D(:,:,j,:)=D;
                end

                %tenBack,tenTar] = obj.TT_RPCA(tensor_4D);
                [tenBack,tenTar] = obj.TR_RPCA(tensor_4D);

                for j = 1: obj.L
                    cur_idx = (i-1) * obj.L + j + obj.start_idx - 1;

                    tarImg = obj.res_patch_ten_mean (tenTar(:,:,j,:), img, obj.patchSize, obj.slideStep);
                    bacImg = obj.res_patch_ten_mean (tenBack(:,:,j,:), img, obj.patchSize, obj.slideStep);

                    T = tarImg / (max(tarImg(:)) + eps);
                    % B = bacImg / (max(bacImg(:)) + eps);
                    tar_res_path = [inpng_res_base_path, num2str(cur_idx,'%04d'), '.png'];
                    imwrite(T, tar_res_path);
                    %obj.result = T;

                end
            end
        end

        function [B, T] = TR_RPCA(obj,D)
            %% initialize parameters
            tol = 1e-7;
            lambda1 = 0;
            lambda2 = 100;
            maxiter = 300;
            H = 2;

            Nway = size(D);
            N = ndims(D);
            L = ceil(N/2);
            sk = zeros(L,1);
            RC = nan(maxiter,2);
            delta = 1./L*ones(1,L);

            epsilon = 1e-3;
            beta = 1e-4;
            gamma = 1e-4;
            beta  = beta * ones(1,L);

            %% initialization
            B = D;
            T = zeros(Nway);
            No = zeros(Nway);
            E = zeros(Nway);
            U = cell(L,1);
            C = cell(L,1);

            for n = 1:L
                U{n}=B;
                C{n}=zeros(Nway);
                order=[n:N 1:n-1];
                lambda1 =lambda1 + H/sqrt(max([prod(Nway(order(1:L))),prod(Nway(order(L+1:N)))]));
                sk(n) = min([prod(Nway(order(1:L))),prod(Nway(order(L+1:N)))]);
            end


            %% ADMM algorithm
            for iter=1:maxiter
                preT = sum(T(:)>0);

                % update U^(n) (auxiliary variables of low-rank part)
                u_s = zeros(Nway);
                c_s = zeros(Nway);

                for n = 1:L
                    order = [n:N 1:n-1];
                    m = permute(B - C{n}/beta(n), order);
                    M = reshape(m,prod(Nway(order(1:L))),[]);
                    [M,sk(n)] = obj.shrink_matrix(M, delta(n)/beta(n), sk(n), epsilon, false);
                    m = reshape(M, Nway(order));
                    U{n} = ipermute(m, order);
                    u_s = u_s + U{n};
                    c_s = c_s + C{n};
                end

                % update B (low-rank part)
                B = (u_s + c_s/gamma+(D - T - No - E/gamma))./(L+1);

                % update T(sparse part)
                T = obj.shrink_vector((D - B - No - E/gamma),lambda1/gamma);

                % update C^(n)
                for n = 1 : L
                    C{n}= C{n} + beta(n) * (U{n} - B);
                end

                % update E
                E = E + gamma * (B + T + No - D);

                %update No
                No = gamma * (D - B - T + E/gamma)/(gamma + 2 * lambda2);

                %update gamma and beta
                gamma = min(1.2 * gamma,1e10);
                beta = min(1.2 * beta,1e10);

                %% Calculate relative error
                currT = sum(T(:) > 0);
                dY = D - B - T - No;
                errList = norm(dY(:)) / norm(D(:));
                disp(['iter ' num2str(iter)  ...
                    ', err=' num2str(errList)...
                    ',|T|0 = ' num2str(currT)]);

                if errList < tol ||(preT>0 && currT>0 && preT == currT)
                    break;
                end

            end

        end

        function [lambda] = weightTC(obj,Nway)
            N = length(Nway);
            lambda = zeros(1,N-1);
            IL = Nway(1);
            for k = 1:N-1
                IR = prod(Nway(k+1:end));
                lambda(k) = min(IL,IR);
                IL = IL*Nway(k+1);
            end
            lambda = lambda/(sum(lambda));
        end

        function [ X ] = SVT(obj, A, tau)

            [U,S,V]=svd(A,'econ');
            s=diag(S);
            idx=find(s>tau,1,'last');
            X=U(:,1:idx)*diag(s(1:idx)-tau)*V(:,1:idx)';
        end

        function [dx] = soft_shrink(obj,v0,tau)
            % s = abs( dx );
            % dx = max(s - Thresh,0).*sign(dx);
            dx=sign(v0).*max(abs(v0)-tau,0);
        end

        function recImg = res_patch_ten_mean(obj, patchTen, img, patchSize, slideStep)

            % 2017-07-31
            % This matlab code implements the RIPT model for infrared target-background
            % separation.
            %
            % Yimian Dai. Questions? yimian.dai@gmail.com
            % Copyright: College of Electronic and Information Engineering,
            %            Nanjing University of Aeronautics and Astronautics

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
                slideStep = 10;
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

        function v1=shrink_vector(obj,v0,tau)
            v1=sign(v0).*max(abs(v0)-tau,0);
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
    end
end