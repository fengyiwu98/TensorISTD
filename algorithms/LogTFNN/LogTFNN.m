% This code is an implement of LogTFNN.
% Ref: X. Kong, C. Yang, S. Cao, C. Li and Z. Peng, "Infrared Small Target Detection via Nonconvex Tensor Fibered Rank Approximation," in IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-21, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: 
% Email:
% Copyright: University of Electronic Science and Technology of China
% Date: 2025/5/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef LogTFNN

    properties
        % para
        patchSize = 40;
        slideStep = 40;

        result;              % tar
    end



    methods
        function obj = process(obj, inImg)

            img = inImg* 256;
            %% constrcut patch tensor of original image
            tenD = obj.gen_patch_ten(img, obj.patchSize, obj.slideStep);
            [n1,n2,n3] = size(tenD);

            %% calculate prior weight map
            %      step 1: calculate two eigenvalues from structure tensor
            [lambda1, lambda2] = obj.structure_tensor_lambda(img, 3);
            %      step 2: calculate corner strength function
            cornerStrength = (((lambda1.*lambda2)./(lambda1 + lambda2)));
            %      step 3: obtain final weight map
            maxValue = (max(lambda1,lambda2));
            priorWeight = mat2gray(cornerStrength .* maxValue);
            %subplot(222),imshow(priorWeight,[]),title('priorWeight image')
            %      step 4: constrcut patch tensor of weight map
            tenW = obj.gen_patch_ten(priorWeight, obj.patchSize, obj.slideStep);

            %% The proposed model
            Nway = size(tenD);
            sigma=1;

            opts=[];
            opts.sigma  = sigma;
            opts.theta  = 0.00001;   % controls \alpha     weight of fibered-rank
            opts.varpi  = 0.15;   % controls \lambda_2  \|T\|_1
            opts.omega  = 10000;   % controls \tau       SVT
            opts.logtol = 0.000001;      % \varepsilon in LogTFNN

            [tenB,tenT,~,~]=obj.my_LogTFNN(tenD,tenW,opts);
            tarImg = obj.res_patch_ten_mean(tenT, img, obj.patchSize, obj.slideStep);
            backImg = obj.res_patch_ten_mean(tenB, img, obj.patchSize, obj.slideStep);
            obj.result = tarImg;

        end

        function patchTen = gen_patch_ten(obj, img, patchSize, slideStep)

            if ~exist('patchSize', 'var')
                patchSize = 50;
            end

            if ~exist('slideStep', 'var')
                slideStep = 10;
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
        end

        function [lambda_1, lambda_2] = structure_tensor_lambda(obj, img, sz)

            G = fspecial('gaussian', [sz sz], 2); % Gaussian kernel
            u = imfilter(img, G, 'symmetric');
            [Gx, Gy] = gradient(u);

            K = fspecial('gaussian', [sz sz], 9); % Gaussian kernel
            J_11 = imfilter(Gx.^2, K, 'symmetric');
            J_12 = imfilter(Gx.*Gy, K, 'symmetric');
            J_21 = J_12;
            J_22 = imfilter(Gy.^2, K, 'symmetric');

            sqrt_delta = sqrt((J_11 - J_22).^2 + 4*J_12.^2);
            lambda_1 = 0.5*(J_11 + J_22 + sqrt_delta);
            lambda_2 = 0.5*(J_11 + J_22 - sqrt_delta);
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

        function [B,T,Z,Out] = my_LogTFNN(obj, X,tenW,opts)

            %% initial value of parameters
            M = obj.rankN(X,0.1);
            % M=100;
            Nway = size(X);

            tol      = 1e-4;
            max_iter = 20;
            rho      = 1.2;

            alpha   = [1,1,opts.theta]/(2+opts.theta);

            lambda0 = 1/sqrt(max(Nway(3),Nway(2))*Nway(1));
            lambda2  = opts.varpi*lambda0;

            tau      = opts.omega*[1,1,1];
            mu       = alpha./tau;
            beta     = 1/mean(tau);

            max_mu   = 1e10*[1,1,1];
            max_beta = 1e10;

            logtol   = opts.logtol;

            %% Initialization
            D=X;
            B = zeros(Nway);
            T = B;
            X1 = B;
            X2 = B;
            X3 = B;
            M1 = X1;
            M2 = X2;
            M3 = X3;
            P = B;
            Z = B;

            V1=ones(Nway);
            V2=ones(Nway);
            V3=ones(Nway);


            Y1=ones(Nway);
            Y2=ones(Nway);
            Y3=ones(Nway);
            Y4=ones(Nway);
            weightTen = ones(Nway);
            Out.Res=[]; Out.PSNR=[];
            for iter = 1 : max_iter
                %% Let
                Lold = B;
                B1 = permute(B,[2,3,1]);  B2 = permute(B,[3,1,2]);  B3 = B;
                m1 = permute(M1,[2,3,1]); m2 = permute(M2,[3,1,2]); m3 = M3;

                %% update X
                tau = alpha./mu;
                X1 = ipermute(obj.ProTlogSum(B1+m1/mu(1),tau(1),logtol),[2,3,1]);
                X2 = ipermute(obj.ProTlogSum(B2+m2/mu(2),tau(2),logtol),[3,1,2]);
                X3 = obj.ProTlogSum(B3+m3/mu(3),tau(3),logtol);

                %% update B
                temp = mu(1)*(X1-M1/mu(1)) + mu(2)*(X2-M2/mu(2)) + mu(3)*(X3-M3/mu(3)) + beta*(D-T+P/beta);
                B = temp/(beta+sum(mu));

                %% update T
                T = obj.prox_l1(D-B+P/beta,weightTen*lambda2/beta);
                weightTen = M./ (abs(T) + 0.01)./tenW;

                %% update Z
                Z=obj.prox_z(Z,Y1,Y2,Y3,Y4,B,rho,V1,V2,V3);

                %% update V1,V2,V3
                dim = size(X);
                tenX = reshape(X, dim);
                dfz1=diff(Z, 1, 1);
                Z1 = zeros(dim);
                Z1(1:end-1,:,:) = dfz1;
                Z1(end,:,:) = tenX(1,:,:) - tenX(end,:,:);

                dfz2=diff(Z, 1, 2);
                Z2 = zeros(dim);
                Z2(:,1:end-1,:) = dfz2;
                Z2(:,end,:) = tenX(:,1,:) - tenX(:,end,:);

                dfz3=diff(Z, 1, 3);
                Z3 = zeros(dim);
                Z3(:,:,1:end-1) = dfz3;
                Z3(:,:,end) = tenX(:,:,1) - tenX(:,:,end);

                dV1=Z1-Y2/rho;
                dV2=Z2-Y3/rho;
                dV3=Z3-Y4/rho;

                V1=obj.prox_l1(dV1, beta/rho);
                V2=obj.prox_l1(dV2, beta/rho);
                V3=obj.prox_l1(dV3, beta/rho);

                %% update Y1,Y2,Y3,Y4
                Y1=Y1+rho*(Z-B);
                Y2=Y2+rho*(V1-Z1);
                Y3=Y3+rho*(V2-Z2);
                Y4=Y4+rho*(V3-Z3);

                %% check the convergence
                dM = D-B-T;
                chg=norm(Lold(:)-B(:))/norm(Lold(:));
                Out.Res = [Out.Res,chg];
                if isfield(opts, 'Xtrue')
                    XT=opts.Xtrue;
                    psnr = PSNR3D(XT * 255, L * 255);
                    Out.PSNR = [Out.PSNR,psnr];
                end

                if mod(iter, 10) == 0
                    if isfield(opts, 'Xtrue')
                        fprintf('3DLogTNN: iter = %d   PSNR= %f   res= %f \n', iter, psnr, chg);
                    else
                        fprintf('3DLogTNN: iter = %d   res= %f \n', iter, chg);
                    end
                end

                if chg < tol
                    break;
                end

                %% update M & P
                P = P + beta*dM;
                M1 = M1 + mu(1)*(B-Z1);
                M2 = M2 + mu(2)*(B-Z2);
                M3 = M3 + mu(3)*(B-Z3);
                beta = min(rho*beta,max_beta);
                mu = min(rho*mu,max_mu);
            end
        end

        function N = rankN(obj, X, ratioN)
            [~,~,n3] = size(X);
            D = obj.Unfold(X,n3,1);
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

        function [X] = ProTlogSum(obj, Y, rho, tol)

            dim = ndims(Y);
            [n1, n2, n3] = size(Y);
            n12 = min(n1, n2);
            Yf = fft(Y, [], dim);
            Uf = zeros(n1, n12, n3);
            Vf = zeros(n2, n12, n3);
            Sf = zeros(n12,n12, n3);

            Yf(isnan(Yf)) = 0;
            Yf(isinf(Yf)) = 0;

            trank = 0;
            endValue = n3/2 + 1;
            for i = 1 : endValue
                [ temp, Sf(:, :, i),Uf(:,:,i),Vf(:,:,i)] = obj.Pro2MlogSum(Yf(:,:,i), rho, tol);
                trank = max(temp, trank);
            end

            for j =n3:-1:endValue+1
                Uf(:,:,j) = conj(Uf(:,:,n3-j+2));
                Vf(:,:,j) = conj(Vf(:,:,n3-j+2));
                Sf(:,:,j) = Sf(:,:,n3-j+2);
            end

            Uf = Uf(:, 1:trank, :);
            Vf = Vf(:, 1:trank, :);
            Sf = Sf(1:trank, 1:trank, :);

            U = ifft(Uf, [], dim);
            S = ifft(Sf, [], dim);
            V = ifft(Vf, [], dim);

            X = obj.tprod( obj.tprod(U,S), obj.tran(V) );
        end

        function [n, SigmaNew, U, V] = Pro2MlogSum(obj, Z, tau, tol)

            [U, Sigma, V] = svd(Z, 'econ');
            Sigma         = diag(Sigma);
            c1            = Sigma-tol;
            c2            = c1.^2-4*(tau-tol*Sigma);
            ind           = find (c2>0);
            n             = length(ind);
            SigmaNew      = zeros(length(Sigma),1) ;
            SigmaNew(1:n) = (c1(1:n)+sqrt(c2(1:n)))/2;
            SigmaNew      = diag(SigmaNew);
        end

        function C = tprod(obj,A,B)

            [n1,~,n3] = size(A);
            l = size(B,2);
            A = fft(A,[],3);
            B = fft(B,[],3);
            C = zeros(n1,l,n3);
            for i = 1 : n3
                C(:,:,i) = A(:,:,i)*B(:,:,i);
            end
            C = real(ifft(C,[],3));
        end

        function Xt = tran(obj,X)

            [n1,n2,n3] = size(X);
            Xt = zeros(n2,n1,n3);
            Xt(:,:,1) = X(:,:,1)';
            if n3 > 1
                for i = 2 : n3
                    Xt(:,:,i) = X(:,:,n3-i+2)';
                end
            end
        end

        function x = prox_l1(obj,b,lambda)

            x = max(0,b-lambda)+min(0,b+lambda);
            x = max(x,0);
        end

        function [X] = prox_z(obj,Z,Y1,Y2,Y3,Y4,L,tho,dV1,dV2,dV3)

            sizeD= size(Z);
            tenX   = reshape(Z, sizeD);
            h = sizeD(1);
            w = sizeD(2);
            d = sizeD(3);
            %%
            Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
            Eny_y   = ( abs(psf2otf([+1, -1], [h,d,w])) ).^2  ;
            Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
            Eny_y   =  permute(Eny_y, [1 3 2]);
            Eny_z   =  permute(Eny_z, [3 1 2]);
            determ  =  Eny_x + Eny_y + Eny_z;

            dfx1=diff(dV1+Y2/tho, 1, 1);
            dfy1=diff(dV2+Y3/tho, 1, 2);
            dfz1=diff(dV3+Y4/tho, 1, 3);

            v1  = zeros(sizeD);
            v2  = zeros(sizeD);
            v3  = zeros(sizeD);
            v1(1:end-1,:,:) = dfx1;
            v1(end,:,:)     = tenX(1,:,:) - tenX(end,:,:);
            v2(:,1:end-1,:) = dfy1;
            v2(:,end,:)     = tenX(:,1,:) - tenX(:,end,:);
            v3(:,:,1:end-1) = dfz1;
            v3(:,:,end)     = tenX(:,:,1) - tenX(:,:,end);

            numer1=L-Y1/tho+v1+v2+v3;
            z  = real( ifftn( fftn(numer1) ./ (determ + 1) ) );
            X  = reshape(z,sizeD);
        end
    end
end