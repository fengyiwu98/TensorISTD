classdef PSTNN

    properties(SetAccess = private)
        patchSize = 10;
        slideStep = 10;
        lambdaL = 0.9;  %tuning0.7
        % flag 是否显示结果
        flag = 1;
    end
    
    properties
        result;
        tarImg;
        backImg;
    end
    
    methods
        
        % inImg must range from 0 to 1
        function obj = process(obj, inImg)
           %% 需要0~255
            inImg = inImg .* 255;
           %% constrcut patch tensor of original image
            tenD = gen_patch_ten(inImg, obj.patchSize, obj.slideStep);
            [n1,n2,n3] = size(tenD);  

           %% calculate prior weight map
            %      step 1: calculate two eigenvalues from structure tensor
            [lambda1, lambda2] = structure_tensor_lambda(inImg, 3);
            %      step 2: calculate corner strength function
            cornerStrength = (((lambda1.*lambda2) ./ (lambda1 + lambda2)));
            %      step 3: obtain final weight map
            maxValue = (max(lambda1, lambda2));
            priorWeight = mat2gray(cornerStrength .* maxValue);
            %      step 4: constrcut patch tensor of weight map
            tenW = gen_patch_ten(priorWeight, obj.patchSize, obj.slideStep);

           %% The proposed model
            lambda = obj.lambdaL / sqrt(max(n1,n2)*n3); 
            [tenB,tenT] = trpca_pstnn(tenD,lambda,tenW); 

           %% recover the target and background image
            obj.tarImg = res_patch_ten_mean(tenT, inImg, obj.patchSize, obj.slideStep);
            obj.backImg = res_patch_ten_mean(tenB, inImg, obj.patchSize, obj.slideStep);
            obj.result = obj.tarImg;
        end
        
    end
end