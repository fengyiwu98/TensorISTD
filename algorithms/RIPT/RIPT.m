% This code is an implement of Reweighted Infrared Patch Tensor.
% Ref: Y. Dai and Y. Wu, "Reweighted Infrared Patch-Tensor Model With Both 
% Nonlocal and Local Priors for Single-Frame Small Target Detection," 
% in IEEE Journal of Selected Topics in Applied Earth Observations and Remote 
% Sensing, vol. 10, no. 8, pp. 3752-3767, Aug. 2017.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Tianfang Zhang
% Email: sparkcarleton@gmail.com
% Copyright:  University of Electronic Science and Technology of China
% Date: 2018/11/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef RIPT
    
    properties
        len = 30;
        step = 10;
        
        loss;
        result;
        
    end
    
    methods
    
        function obj = process(obj, inImg)
% tic
            [m, n] = size(inImg);

            % lambda = 0.7 / sqrt(min(m, n));
            lambda = 1 / (length(1:obj.step:m-obj.len)*length(1:obj.step:n-obj.len));
            h = 1;

            % Calculate the local structure weight
            W_LS = obj.LocalStructure(inImg, 3, h);
%             figure('numbertitle', 'off', 'name', 'Local Structure Weight');
%             imshow(W_LS ./ max(W_LS(:)), []);

            % Convert the image and local structure weight to tensor
            tenF = zeros(obj.len, obj.len, length(1:obj.step:m-obj.len)*length(1:obj.step:n-obj.len));
            tenWLS = zeros(size(tenF));
            counter = 1;
            for i = 1:obj.step:m-obj.len
                for j = 1:obj.step:n-obj.len
                    tenF(:, :, counter) = inImg(i:i+obj.len-1, j:j+obj.len-1);
                    tenWLS(:, :, counter) = W_LS(i:i+obj.len-1, j:j+obj.len-1);
                    counter = counter + 1;
                end
            end

            % Show the singulars in all mode
%             ShowSingular(tenF);

            % Solve the robust PCA problem via ADMM
            [tenT, tenB, obj.loss] = obj.func_ADMM(tenF, tenWLS, lambda);

            % Convert the tensor to image
            reconsT = zeros(m, n, 100);
            reconsB = zeros(m, n, 100);
            countMatrix = zeros(m, n);
            index = 1;
            for i = 1:obj.step:m-obj.len
                for j = 1:obj.step:n-obj.len
                    % Count the time each pixel used and record the value
                    repatch_T = tenT(:, :, index);
                    repatch_B = tenB(:, :, index);
                    countMatrix(i:i+obj.len-1, j:j+obj.len-1) = countMatrix(i:i+obj.len-1, j:j+obj.len-1) + 1;

                    % Record the value of each pixel
                    for ii = i:i+obj.len-1
                        for jj = j:j+obj.len-1
                            reconsT(ii, jj, countMatrix(ii,jj)) = repatch_T(ii-i+1, jj-j+1);
                            reconsB(ii, jj, countMatrix(ii,jj)) = repatch_B(ii-i+1, jj-j+1);
                        end
                    end

                    index = index + 1;
                end
            end

            rstT = zeros(m, n);
            rstB = zeros(m, n);

            for i = 1:m
                for j = 1:n
                    % Mean
                    if countMatrix(i ,j) > 0
                        vectorT = reconsT(i, j, 1:countMatrix(i,j));
                        vectorB = reconsB(i, j, 1:countMatrix(i,j));

                        rstT(i, j) = median(vectorT);
                        rstB(i, j) = median(vectorB);
                    end
                end
            end
            rstT = rstT .* (rstT > 0);
            % rstT = rstT / max(rstT(:));
            rstB = rstB .* (rstB > 0);
% toc
            obj.result = rstT;
            
%             % Show the result
%             figure('numbertitle', 'off', 'name', 'The result');
%             subplot(121),imshow(rstT, []),title('Target');
%             subplot(122),imshow(rstB, []),title('Background');
%             [X, Y] = meshgrid(1:n, 1:m);
%             figure,mesh(X, Y, rstT);
        end
        
        function [tenT, tenB, loss] = func_ADMM(obj, tenF, WLS, lambda)

            % Parameter Initialization
            tenT_k = zeros(size(tenF));
            tenYs_k = cell(1, 3);   tenYs_kp1 = cell(1, 3);
            tenBs_k = cell(1, 3);   tenBs_kp1 = cell(1, 3);
            for i = 1:3
                tenBs_k{i} = tenF;
                tenYs_k{i} = zeros(size(tenF));
            end
            Weight_k = ones(size(tenF)) .* WLS;
            mu_k = 5 * std(tenF(:));

            normF = norm(tenF(:), 'fro');
            rho = 1.5;
            iterNum = 0;
            Converged = 0;
            loss = [];


            % Iteration steps
            while ~Converged
                % Update background tensors Bi
                for i = 1:3
                    mat_tmp = obj.ten2mat(tenF + mu_k*tenYs_k{i} - tenT_k, i);
                    mat_tmp = obj.singularValueShrinkage(mat_tmp, mu_k);
                    tenBs_kp1{i} = obj.mat2ten(mat_tmp, size(tenF), i);
                end

                % Update target tensors T
                ten_tmp_mean = zeros(size(tenF));
                for i = 1:3
                    ten_tmp_mean = ten_tmp_mean + tenF + mu_k*tenYs_k{i} - tenBs_kp1{i};
                end
                ten_tmp_mean = ten_tmp_mean / 3;
                tenT_kp1 = obj.softThreshold(ten_tmp_mean .* (ten_tmp_mean>0), mu_k * lambda * Weight_k / 3);

                % Update lagrange multiplier Yi
                for i = 1:3
                    tenYs_kp1{i} = tenYs_k{i} + (tenF - tenBs_kp1{i} - tenT_kp1)/mu_k;
                end

                % Update other parameters
                Weight_kp1 = WLS ./ (abs(tenT_kp1) + 0.01);
                mu_kp1 = mu_k / rho;
                iterNum = iterNum + 1;

                % Judge converged
                tenB_check = zeros(size(tenF));
                for i = 1:3
                    tenB_check = tenB_check + tenBs_kp1{i};
                end
                tenB_check = tenB_check / 3;
                stopCriterion = norm(tenF(:) - tenB_check(:) - tenT_kp1(:), 'fro') / normF;
                loss(iterNum) = stopCriterion;
                nonzeroNum_k = sum(tenT_k(:) ~= 0);
                nonzeroNum_kp1 = sum(tenT_kp1(:) ~= 0);

                if stopCriterion < 1e-7 || nonzeroNum_k == nonzeroNum_kp1
                    Converged = 1;
                end


                % Display the status
%                 disp(['Iteration Number: ', num2str(iterNum), '    loss: ', num2str(stopCriterion)]);
%                 disp(['nonzeroNum_k: ', num2str(nonzeroNum_k)]);
%                 disp(' ');

                % Assignment to continue loop
                tenT_k = tenT_kp1;
                Weight_k = Weight_kp1;
                mu_k = mu_kp1;
                for i = 1:3
                    tenBs_k{i} = tenBs_kp1{i};
                    tenYs_k{i} = tenYs_kp1{i};
                end

            end

            % Assignment
            tenB_tmp = zeros(size(tenF));
            for i = 1:3
                tenB_tmp = tenB_tmp + tenBs_k{i};
            end
            tenB = tenB_tmp / 3;
            tenT = tenT_k;

        end
        
        function W_LS = LocalStructure(obj, inImg, filter_size, h)

            [m, n] = size(inImg);
            gaussian_filter = fspecial('gaussian', [filter_size, filter_size], 1.2);
            u = imfilter(inImg, gaussian_filter, 'replicate');

            [Gx, Gy] = gradient(u);
            J = zeros(m, n, 2, 2);

            % Define another filter
            K = fspecial('gaussian', [filter_size, filter_size], 2);
            J(:, :, 1, 1) = imfilter(Gx.*Gx, K, 'replicate');
            J(:, :, 1, 2) = imfilter(Gx.*Gy, K, 'replicate');
            J(:, :, 2, 1) = J(:, :, 1, 2);
            J(:, :, 2, 2) = imfilter(Gy.*Gy, K, 'replicate');

            tmp = sqrt((J(:,:,2,2)-J(:,:,1,1)).^2 + 4*J(:,:,1,2).^2);
            lambda_1 = J(:,:,1,1) + J(:,:,2,2) + tmp;
            lambda_2 = J(:,:,1,1) + J(:,:,2,2) - tmp;

            dmax = max(max(lambda_1 - lambda_2));
            dmin = min(min(lambda_1 - lambda_2));

            W_LS = exp(h * (mat2gray(lambda_1 - lambda_2)));
    
        end
        
        function ten = mat2ten(obj, mat, size_arr, i)

            % 2017-07-31
            % This matlab code implements the RIPT model for infrared target-background 
            % separation.
            %
            % Yimian Dai. Questions? yimian.dai@gmail.com
            % Copyright: College of Electronic and Information Engineering, 
            %            Nanjing University of Aeronautics and Astronautics  

            perm_size = size_arr;
            perm_size(i) = [];
            perm_size = [size_arr(i), perm_size];
            perm_ten = reshape(mat, perm_size);
            N = length(size_arr);
            order = 2:N;
            if i == N
                order = [order, 1];
            elseif i == 1
                order = [1, order];
            else       
                order = [order(1:i-1), 1, order(i:end)];    
            end
            ten = permute(perm_ten, order);
        end
        
        function output = singularValueShrinkage(obj, inImg, epsilon)

            [U, S, V] = svd(inImg, 'econ');
            output = U*obj.softThreshold(S, epsilon)*V';

        end
        
        function output = softThreshold(obj, inImg, epsilon)
            output = zeros(size(inImg));
            posImg = inImg .* (inImg>epsilon);
            negImg = inImg .* (inImg<-epsilon);
            output = output + (posImg-epsilon).*(inImg>epsilon) + (negImg+epsilon) .* (inImg<-epsilon);
        end
        
        function mat = ten2mat(obj, ten, i)

            % 2017-07-31
            % This matlab code implements the RIPT model for infrared target-background 
            % separation.
            %
            % Yimian Dai. Questions? yimian.dai@gmail.com
            % Copyright: College of Electronic and Information Engineering, 
            %            Nanjing University of Aeronautics and Astronautics  

            N = ndims(ten);
            order = 1:N;
            order(i) = [];
            order = [i, order];
            perm_ten = permute(ten, order);
            ten_size = size(perm_ten);
            rows = ten_size(1);
            cols = prod(ten_size(2:end));
            mat = reshape(perm_ten, [rows, cols]);
        end
    
    end
    
end