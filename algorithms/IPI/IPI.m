classdef IPI

    properties
        opt;
        loss;
        result;
    end
    
    methods
        
        function obj = IPI
            obj.opt = struct();
            obj.opt.dw = 60;%30
            obj.opt.dh = 60;
            obj.opt.x_step = 25;
            obj.opt.y_step = 25;%10
        end
        
        function [A_hat, E_hat, loss] = winRPCA_median(obj, I, opt)

            if size(I, 3) == 3
                I = rgb2gray(I);
            end
            [m n] = size(I);
            A_hat = zeros(m,n);
            E_hat = zeros(m ,n);
            C = zeros(m, n);
            dw = opt.dw;
            dh = opt.dh;
            x_step = opt.x_step;
            y_step = opt.y_step;

            I1 = I ; %imfilter(I, fspecial('gaussian', 5));
            D = [];% 块图像 patch-image
            for i = 1:y_step:m-dh+1
                for j = 1:x_step:n-dw+1
                    temp = I1(i:i+dh-1, j:j+dw-1);
                    D = [D, temp(:)];
                end
            end
            D = mat2gray(D);
            % D = D - mean2(D);
            [m1 n1] = size(D);
            lambda = 1 / sqrt(max(m, n));
            
            % 优化问题求解
            [A1, E1, loss] = obj.APG_IR(D, lambda);
            
            %重建恢复
            AA = zeros(m, n, 100);
            EE = zeros(m, n, 100);
            temp = zeros(dh, dw); % 定义成30*30大小的矩阵
            temp1 = zeros(dh, dw);
            index = 0;
            for i = 1:y_step:m-dh+1 % 每个块30*30的左上角 即对块进行循环访问
                for j = 1:x_step:n-dw+1
                    index = 1+index;
                    temp(:) = A1(:, index); % 每个块,即分解优化后的A1和E1的一列
                    temp1(:) = E1(:, index); % matlab默认的赋值是先列后行,将1*900的列向量依次赋值到30*30矩形中
                    C(i:i+dh-1, j:j+dw-1) = C(i:i+dh-1, j:j+dw-1) + 1; % C记录了原图某个像素点是属于那个块 即块的编号 从1开始
                    for ii = i:i+dh-1 % 每个块内的像素处理 A1和E1每列进行reshape原来块的形状30*30
                        for jj = j:j+dw-1  % ii-i+1 计算的是当前处理的坐标在块patch中的相对位置
                            AA(ii,jj, C(ii,jj)) = temp(ii-i+1, jj-j+1);
                            EE(ii,jj, C(ii,jj)) = temp1(ii-i+1, jj-j+1);
                        end
                    end
                end
            end
            %     C(find(C==0)) = 1000;
            for i=1:m
                for j=1:n
                    if C(i,j) > 0
                        A_hat(i,j) = median(AA(i,j,1:C(i,j))); % 
                        E_hat(i,j) = median(EE(i,j,1:C(i,j)));
                    end
                end
            end
        end
        
        function [A_hat,E_hat, loss] = APG_IR(obj, D, lambda, maxIter, tol, ...
            lineSearchFlag, continuationFlag, eta, mu, outputFileName)
            % This code is an implementation of Proximal Gradient Algorithm for small target
            % detecton in our published paper: Chenqiang Gao, Deyu Meng, Yi Yang, et al., "Infrared Patch-Image Model for Small Target Detection in a Single Image, " Image Processing, IEEE Transactions on, vol. 22, no. 12, pp. 4996-5009, 2013.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If you use this code in your publications, please cite:
            % Chenqiang Gao, Deyu Meng, Yi Yang, et al., "Infrared Patch-Image Model for Small Target Detection in a Single Image, " Image Processing, IEEE Transactions on, vol. 22, no. 12, pp. 4996-5009, 2013. 
            % @article{Gao2013,
            %    author = {Gao, Chenqiang and Meng, Deyu and Yang, Yi and Wang, Yongtao and Zhou, Xiaofang and Hauptmann, Alex},
            %    title = {Infrared Patch-Image Model for Small Target Detection in a Single Image},
            %    journal = {Image Processing, IEEE Transactions on},
            %    volume = {22},
            %    number = {12},
            %    pages = {4996-5009},
            %    year = {2013}
            % }
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This code is implemented based on the original code of Accelerated
            % Proximal Gradient algorithm on the websie:
            % "http://perception.csl.illinois.edu/matrix-rank/sample_code.html". So
            % please also cite the following reference:
            % "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
            % "An Accelerated Proximal Gradient Algorithm for Nuclear Norm Regularized Least Squares problems", K.-C. Toh and S. Yun, preprint 2009.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Paremeter Explanation
            %
            % D - m x n matrix of observations/data (required input)
            % lambda - weight on sparse error term in the cost function (required input)
            %
            % tol - tolerance for stopping criterion.
            %     - DEFAULT 1e-7 if omitted or -1.
            % maxIter - maximum number of iterations
            %         - DEFAULT 10000, if omitted or -1.
            % lineSearchFlag - 1 if line search is to be done every iteration
            %                - DEFAULT 0, if omitted or -1.
            % continuationFlag - 1 if a continuation is to be done on the parameter mu
            %                  - DEFAULT 1, if omitted or -1.
            % eta - line search parameter, should be in (0,1)
            %     - ignored if lineSearchFlag is 0.
            %     - DEFAULT 0.9, if omitted or -1.
            % mu - relaxation parameter
            %    - ignored if continuationFlag is 1.
            %    - DEFAULT 1e-3, if omitted or -1.
            % outputFileName - Details of each iteration are dumped here, if provided.
            %
            % [A_hat, E_hat] - estimates for the low-rank part and error part, respectively
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If you have any questions, please contact:
            % Author: Chenqiang Gao 
            % Email: gaochenqiang@gmail.com or gaocq@cqupt.edu.cn
            % Copyright:  Chongqing University of Posts and Telecommunications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if nargin < 2
                error('Too few arguments') ;
            end

            if nargin < 4
                tol = 1e-6 ;
            elseif tol == -1
                tol = 1e-6 ;
            end

            if nargin <= 3
                maxIter = 1000 ;
            elseif maxIter == -1
                maxIter = 1000 ;
            end

            if nargin < 6
                continuationFlag = 1 ;
            elseif continuationFlag == -1 ;
                continuationFlag = 1 ;
            end

            if ~continuationFlag
                if nargin < 8
                    mu = 1e-3 ;
                elseif mu == -1
                    mu = 1e-3 ;
                end
            end
            % mu = 1e6 ;

            DISPLAY_EVERY = 20 ;
            DISPLAY = 0 ;
            
            loss = [];

            % Initializing optimization variables

            [m, n] = size(D) ;
            t_k = 1 ; % t^k
            t_km1 = 1 ; % t^{k-1}

            tau_0 = 2 ; % square of Lipschitz constant for the RPCA problem

            X_km1_A = zeros(m,n) ; X_km1_E = zeros(m,n) ; % X^{k-1} = (A^{k-1},E^{k-1})
            X_k_A = zeros(m,n) ; X_k_E = zeros(m,n) ; % X^{k} = (A^{k},E^{k})


            if continuationFlag
                mu_0 = norm(D) ;
                mu_k = 0.99*mu_0 ;
                mu_bar = 1e-9 * mu_0 ;
            else
                mu_k = mu ;
            end
            [U S V] = svd(D, 'econ');
            s = diag(S);
            mu_k = s(2);
            mu_bar = 0.005 * s(4);

            tau_k = tau_0 ;

            converged = 0 ;
            numIter = 0 ;

            NOChange_counter = 0;
            pre_rank = 0;
            pre_cardE = 0;
            % tol = 1e-6 ;

            % Start main loop
            while ~converged

                Y_k_A = X_k_A + ((t_km1 - 1)/t_k)*(X_k_A-X_km1_A) ;
                Y_k_E = X_k_E + ((t_km1 - 1)/t_k)*(X_k_E-X_km1_E) ;

                G_k_A = Y_k_A - (1/tau_k)*(Y_k_A+Y_k_E-D) ;
                G_k_E = Y_k_E - (1/tau_k)*(Y_k_A+Y_k_E-D) ;

                [U S V] = svd(G_k_A, 'econ');
                diagS = diag(S) ;
            %      temp = max(diagS(150), mu_k/tau_k);

            %     X_kp1_A = U * diag(pos(diagS - temp)) * V';
                X_kp1_A = U * diag(obj.pos(diagS- mu_k/tau_k)) * V';

            %     X_kp1_E = sign(G_k_E) .* pos( abs(G_k_E) - lambda* temp );
                X_kp1_E = sign(G_k_E) .* obj.pos( abs(G_k_E) - lambda* mu_k/tau_k );

            %   rankA  = sum(diagS>temp);
                rankA  = sum(diagS>mu_k/tau_k);

                cardE = sum(sum(double(abs(X_kp1_E)>0)));

                t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;

                temp = X_kp1_A + X_kp1_E - Y_k_A - Y_k_E ;
                S_kp1_A = tau_k*(Y_k_A-X_kp1_A) + temp ;
                S_kp1_E = tau_k*(Y_k_E-X_kp1_E) + temp ;


                stoppingCriterion = norm([S_kp1_A,S_kp1_E],'fro') / (tau_k * max(1, norm([X_kp1_A, X_kp1_E],'fro'))) ;

                if stoppingCriterion <= tol
                    converged = 1 ;
                end

                if continuationFlag
                    mu_k = max(0.9*mu_k, mu_bar) ;
                end

                t_km1 = t_k ;
                t_k = t_kp1 ;
                X_km1_A = X_k_A ; X_km1_E = X_k_E ;
                X_k_A = X_kp1_A ; X_k_E = X_kp1_E ;

                numIter = numIter + 1 ;
                loss(numIter) = stoppingCriterion;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % The iteration process can be finished if the rank of A keeps the same
                % many times;
                if pre_rank == rankA
                    NOChange_counter = NOChange_counter +1;
                    if NOChange_counter > 10 && abs(cardE-pre_cardE) < 20
                        converged = 1 ;
                    end
                else
                    NOChange_counter = 0;
                    pre_cardE = cardE;
                end
                pre_rank = rankA;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % In practice, the APG algorithm, sometimes, cannot get  
                % a strictly low-rank matrix A_hat after iteration process. Many
                % singular valus of the obtained matrix A_hat, however, are extremely small. This can be considered 
                % to be low-rank to a certain extent. Experimental results show that the final recoverd 
                % backgournd image and target image are good.
                % Alternatively, we can make the rank of matrix A_hat lower using the following truncation. 
                % This trick can make the APG algorithm faster and the performance of our algorithm is still satisfied. 
                % Here we set the truncated threshold as 0.3, while it can be adaptively set based on your actual scenarios.
                if rankA > 0.3 * min([m n]) 
                    converged = 1 ;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if DISPLAY && mod(numIter,DISPLAY_EVERY) == 0
                    disp(['Iteration ' num2str(numIter) '  rank(A) ' num2str(rankA) ...
                        ' ||E||_0 ' num2str(cardE)])
                end

                if nargin > 8
                    fprintf(fid, '%s\n', ['Iteration ' num2str(numIter) '  rank(A)  ' num2str(rankA) ...
                        '  ||E||_0  ' num2str(cardE) '  Stopping Criterion   ' ...
                        num2str(stoppingCriterion)]) ;
                end

                if ~converged && numIter >= maxIter
                    disp('Maximum iterations reached') ;
                    converged = 1 ;
                end

            end

            A_hat = X_k_A ;
            E_hat = X_k_E ;
        end
        
        
        function P = pos(obj, A)

            P = A .* double( A > 0 );
            
        end
        
        function obj = process(obj, inImg)
%         tic
            [A, E, loss] = obj.winRPCA_median(inImg, obj.opt);
            obj.result = E .* (E>0);%舍弃负数
            obj.loss = loss;
%         toc
%         result = obj.result;
%          % Show the result
%             figure
%             imshow(result, []),title('Target');
        end
    end
    
    
end