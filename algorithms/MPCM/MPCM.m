 % This code is an implement of MPCM.
% Ref: Wei Y, You X, Li H. Multiscale patch-based contrast measure for small infrared target detection[J]. Pattern Recognition, 2016, 58(C):216-226.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Tianfang Zhang
% Email: sparkcarleton@gmail                                    .com
% Copyright:  University of Electronic Science and Technology of China
% Date: 2018/10/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

classdef MPCM

    properties
        radius = 1:9;
        result;
    end
    
    methods
        
        function obj = process(obj, inImg)
            
            [m, n] = size(inImg);
            
            in_diameters = 2*obj.radius + 1;
            ex_diameters = 3 * in_diameters;

            rstImgs = zeros(m, n, length(obj.radius));
            for k = 1:length(obj.radius)    % loop each scale
                % mean filtering
                r = obj.radius(k) + 1;
                d = in_diameters(k);
                meanImg = imfilter(inImg, ones(d,d)/d^2, 'replicate');

                % construct filters

                filter1 = zeros(ex_diameters(k), ex_diameters(k));
                filter1(r, r) = -1;
                filter1(r+d, r+d) = 1;

                filter2 = zeros(ex_diameters(k), ex_diameters(k));
                filter2(r, r+d) = -1;
                filter2(r+d, r+d) = 1;

                filter3 = zeros(ex_diameters(k), ex_diameters(k));
                filter3(r, r+2*d) = -1;
                filter3(r+d, r+d) = 1;

                filter4 = zeros(ex_diameters(k), ex_diameters(k));
                filter4(r+d, r+2*d) = -1;
                filter4(r+d, r+d) = 1;

                filter5 = zeros(ex_diameters(k), ex_diameters(k));
                filter5(r+2*d, r+2*d) = -1;
                filter5(r+d, r+d) = 1;

                filter6 = zeros(ex_diameters(k), ex_diameters(k));
                filter6(r+2*d, r+d) = -1;
                filter6(r+d, r+d) = 1;

                filter7 = zeros(ex_diameters(k), ex_diameters(k));
                filter7(r+2*d, r) = -1;
                filter7(r+d, r+d) = 1;

                filter8 = zeros(ex_diameters(k), ex_diameters(k));
                filter8(r+d, r) = -1;
                filter8(r+d, r+d) = 1;

            %     figure;
            %     subplot(241),imshow(filter1, []),title('filter 1');
            %     subplot(242),imshow(filter2, []),title('filter 2');
            %     subplot(243),imshow(filter3, []),title('filter 3');
            %     subplot(244),imshow(filter4, []),title('filter 4');
            %     subplot(245),imshow(filter5, []),title('filter 5');
            %     subplot(246),imshow(filter6, []),title('filter 6');
            %     subplot(247),imshow(filter7, []),title('filter 7');
            %     subplot(248),imshow(filter8, []),title('filter 8');

                % filtering
                fs = zeros(m, n, 8);
                fs(:, :, 1) = imfilter(meanImg, filter1, 'replicate');
                fs(:, :, 2) = imfilter(meanImg, filter2, 'replicate');
                fs(:, :, 3) = imfilter(meanImg, filter3, 'replicate');
                fs(:, :, 4) = imfilter(meanImg, filter4, 'replicate');
                fs(:, :, 5) = imfilter(meanImg, filter5, 'replicate');
                fs(:, :, 6) = imfilter(meanImg, filter6, 'replicate');
                fs(:, :, 7) = imfilter(meanImg, filter7, 'replicate');
                fs(:, :, 8) = imfilter(meanImg, filter8, 'replicate');

                cs = zeros(m, n, 4);
                cs(:, :, 1) = fs(:, :, 1) .* fs(:, :, 5);
                cs(:, :, 2) = fs(:, :, 2) .* fs(:, :, 6);
                cs(:, :, 3) = fs(:, :, 3) .* fs(:, :, 7);
                cs(:, :, 4) = fs(:, :, 4) .* fs(:, :, 8);

                rstImgs(:, :, k) = min(cs, [], 3); 
            end

            obj.result = max(rstImgs, [], 3);
            obj.result = obj.result .* (obj.result>0);
        end
        
    end

end