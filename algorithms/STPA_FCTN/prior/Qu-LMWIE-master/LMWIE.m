function [out] = LMWIE(img, adjmat, ck)
%LMWIE 
% Matlab realizaion of 'Local Mutation Weighted information Entropy' from
% “Novel detection method for infrared small targets using weighted information entropy,”
%  by X.Qu in 2012.
% In this realization, an improved local window with flexible shape and
% pixel filter criteria is developed based on origin design.
%
% 'adjmat' is a boolean matrix enums pixels on local grid with non-zero 
% value within neighborhoods, which brings flexibility on window size and
% shape.

% convert input frame into double format
if(ischar(img))
    img = imread(img);
end
[imgR ,imgC, dimension]=size(img);
if(dimension>2) 
    img = rgb2gray(img);
    dimension = 1;
end

% LMWIE is suggested to process input image in uint8 format.
%img = double(img);

% parse adjacency matrix into local offset list for slide-sampling
[lwR, lwC] = size(adjmat);
cent_y = fix((lwR+1)/2);
cent_x = fix((lwC+1)/2);

% sampling coord offset from window centre
% 2*N matrix stores y and x offset in row 1 and 2.
smp_ofs_list = zeros(2, lwR*lwC);   
smp_ofs_size = 0;
for y = 1:lwR
    for x = 1:lwC
        % select non-zero position as valid sampling offset 
        if(adjmat(y,x))
            smp_ofs_size = smp_ofs_size + 1;
            smp_ofs_list(1, smp_ofs_size) = y - cent_y;
            smp_ofs_list(2, smp_ofs_size) = x - cent_x;
        end
    end
end
smp_ofs_list = smp_ofs_list(:,1:smp_ofs_size);

% slide adjacency matrix across input image with iterating Local Entropy.
out = zeros(imgR, imgC);
for y = 1:imgR
    for x = 1:imgC
        % sampled pixel from local adjacency of ongoing pixel
        smp_qty = 0;
        tmp_pix_smp = zeros(3, smp_ofs_size); 
        tmp_adj_gray = zeros(1, 256); % gray scales of local neighborhood
        for i = 1:smp_ofs_size
            smp_y = y + smp_ofs_list(1, i);
            smp_x = x + smp_ofs_list(2, i);
            if (smp_y>0 && smp_y<=imgR) && (smp_x>0 && smp_x<=imgC)
                smp_qty = smp_qty + 1;
                pix_gs = img(smp_y,smp_x);
                tmp_pix_smp(1, smp_qty) = smp_y;
                tmp_pix_smp(2, smp_qty) = smp_x;
                tmp_pix_smp(3, smp_qty) = pix_gs;
                gs_idx = pix_gs + 1;
                tmp_adj_gray(gs_idx) = tmp_adj_gray(gs_idx) + 1;
            end
        end
        tmp_pix_smp = tmp_pix_smp(:,1:smp_qty);
        
        % calculate Enhanced pixel of current position via Local Weighted
        % Entropy (See details in 3.1 & 3.2 of the paper).
        [out(y,x), tmp_pix_smp, ~] = pixSeq_LWE(img(y,x), tmp_pix_smp, tmp_adj_gray);
    end
end

bw = out_bw(out, ck);

% visualization
figure(1);
subplot(1,3,1);
imshow(img,[]);
subplot(1,3,2);
imshow(out,[]);
subplot(1,3,3);
imshow(bw,[]);

end

