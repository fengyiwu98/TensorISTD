function [out] = HOG_Saliency(img)
%UNTITLED 此处显示有关此函数的摘要
% A MATLAB realization on 'Infrared small target detection based on saliency
% and gradients difference measure', published on 'Optical and Quantum
% Electronics' 2020(52):151, DOI of which is
% https://doi.org/10.1007/s11082-020-2197-x
% Build 200707, by IraDei.

% convert input frame into double format
if(ischar(img))
    img = imread(img);
end
[imgR ,imgC, dimension]=size(img);
if(dimension>2) 
    img = rgb2gray(img);
    dimension = 1;
end

% Fourrier-transform input image into frequency salient map
img=double(img);

frr = fft2(img);
frr_sft = fftshift(frr);
f_mean_size = 3;    % mean filter size in frequency domain
A = abs(frr);   % frequency Amplitude
L = log(A); % log of frequency Amplitude
P = log(angle(frr)*180/pi); % phase specturm of Fourrier trans
R = L - imfilter(L, fspecial('average', [f_mean_size, f_mean_size]));  % Frequency Saliency map

% Reconstruct spatial saliency with inverse conversion based on frequency
% map and phase spectrum. 
% origin eq.6 given in 3.1 WRONGLY expressed spectral reconstruction with
% merging frequency salency into phase exponent term.
% spatial_saliency = ifft2(exp(R + 1i*P));
spatial_saliency = ifft2(abs(R).*exp(1i*angle(frr)));

sg_size = 3;    % spatial gaussian filter size
ssigma = 2.5;   % sigma of spatial gaussian filter
% square amplitude in reconstructed saliency map to highlight potential targets.
spatial_saliency_sqr = spatial_saliency.*spatial_saliency;
Sa = imfilter(spatial_saliency_sqr, fspecial('gaussian', sg_size, ssigma));

% HoG saliency extraction
% Create the operators for computing image derivative at every pixel.
hx = [-1,0,1];
hy = hx';
% Calculate whole image gradient so we just need to cite data here in later
% HoG sliding sample phase.
% Compute the derivative in the x and y direction for every pixel.
dx = filter2(hx, double(img));
dy = filter2(hy, double(img));
% Convert the gradient vectors to polar coordinates (angle and magnitude).
angles = atan2(dy, dx);
magnit = ((dy.^2) + (dx.^2)).^.5;

HoG_Saliency = zeros(imgR, imgC);
% HoG difference 
% Init HoG Sampling structure slides across the input image
hog = setHoGParams(9, 1, 3); 
smp_npix = hog.winR*hog.winC;
% slide sampling structure with non-padding manner
for y = 1:imgR
    for x = 1:imgC
        % deploy sample kernel and enum pixels within sample neighbor into
        % temporary sampling image.
        
        % writing pointer on temp sample image 
        tmp_smp_ymin = 0;
        tmp_smp_xmin = 0;
        tmp_smp_xptr = NaN;
        tmp_smp_yptr = NaN;
        
        % Init containers for local sample cell 
        smp_cellFlag = zeros(1, hog.ncell);
        smp_cellAngles = cell(hog.numVertCells, hog.numHorizCells);
        smp_cellMagnitudes = cell(hog.numVertCells, hog.numHorizCells);
        % Commence local sampling with coord verification
        for ofs_y = 1:hog.winR
            smp_y = y + hog.smpofs(ofs_y, 1, 1);
             % calculate and verify row potision of sampled pixel on image plane
            if (smp_y>0 && smp_y<=imgR)
                if(tmp_smp_ymin==0)
                    tmp_smp_ymin = ofs_y;
                    tmp_smp_yptr = tmp_smp_ymin;
                end
                % verify col position 
                for ofs_x = 1:hog.winC
                    smp_x = x + hog.smpofs(ofs_y,ofs_x, 2);
                    if(smp_x>0 && smp_x<=imgC)
                        if(tmp_smp_xmin==0)
                            tmp_smp_xmin = ofs_x;
                            tmp_smp_xptr = tmp_smp_xmin;
                        end
                        
                        % allocate cell index and cell position to local sample image 
                        smp_cell_idx = hog.cellmap(tmp_smp_yptr, tmp_smp_xptr, 1);
                        tmp_cell_posy = hog.cellmap(tmp_smp_yptr, tmp_smp_xptr, 2);
                        tmp_cell_posx = hog.cellmap(tmp_smp_yptr, tmp_smp_xptr, 3);
                        % cell position in block matrix
                        smp_cblk_posy = hog.cidx_pos(smp_cell_idx, 1);
                        smp_cblk_posx = hog.cidx_pos(smp_cell_idx, 2);
                        % enum pixel gradient data into from global map
                        if(~smp_cellFlag(smp_cell_idx))
                            % setup cell container for new cell enumed
                            smp_cellFlag(smp_cell_idx) = 1;
                            smp_cellAngles{smp_cblk_posy, smp_cblk_posx} = zeros(hog.cellSize);
                            smp_cellMagnitudes{smp_cblk_posy, smp_cblk_posx} = zeros(hog.cellSize);
                        end
                        
                        % write sample cell index map
                        smp_cellAngles{smp_cblk_posy, smp_cblk_posx}(tmp_cell_posy, tmp_cell_posx) = angles(smp_y, smp_x);
                        smp_cellMagnitudes{smp_cblk_posy, smp_cblk_posx}(tmp_cell_posy, tmp_cell_posx) = magnit(smp_y, smp_x);
                        tmp_smp_xptr = tmp_smp_xptr + 1;
                    elseif(smp_x>imgC)
                        break;
                    end
                end
                
                % Update temporary sampling inwriting pointers
                tmp_smp_yptr = tmp_smp_yptr + 1;
                tmp_smp_xptr = tmp_smp_xmin;
            elseif(smp_y>imgR)
                break;
            end
        end
        
        % calculate HoG of each member cell of whole sample window, thus we
        % could avoid redundant computations at block level with
        % normalization on member cell combination according to block list.
        
        % Init temp HoG container for 'hog.ncell' member cells in sample
        % window.
        tmp_Hist = zeros(hog.numVertCells, hog.numHorizCells, hog.numBins);
        % calculate HoG of each member cell in current sample window
        for i = 1:hog.numVertCells
            for j = 1:hog.numHorizCells
                cidx = (i-1)*hog.numHorizCells + j;
                if(smp_cellFlag(cidx))
                    tmp_Hist(i, j, :) = getHistogram(smp_cellMagnitudes{i, j}, smp_cellAngles{i, j}, hog.numBins);
                end
            end
        end
        % Contrast block in HoG-LCM MUST ALWAYS contain the central target cell T,
        % thus cell-rank of single block ALWAYS equals to halved rank 
        % , hog.rkP, of the cell matrix. This makes total block quantity 
        % ALWAYS equals to cell quantity, hog.nP, in single contrast block.
        tmp_blk_Hist = cell(1, hog.nP);
        tmp_blk_HistN = cell(1, hog.nP);
        tmp_blk_HistN_Bin = zeros(hog.nP, hog.numBins);
        for blk_idx = 1:hog.nP
            % normalize HoG of all member cells in each contrast block
            for j = 1:hog.nP
                % enum HoG of j-th non-empty cell belongs to current contrast block
                cell_idx = hog.block_cell_list(blk_idx, j);
                if(smp_cellFlag(cell_idx))
                    smp_cblk_posy = hog.cidx_pos(cell_idx, 1);
                    smp_cblk_posx = hog.cidx_pos(cell_idx, 2);
                    tmp_blk_Hist{blk_idx} = [tmp_blk_Hist{blk_idx}, tmp_Hist(smp_cblk_posy, smp_cblk_posx, :)];
                end  
            end
            % Put all the histogram values into a single vector (nevermind the 
            % order), and compute the magnitude.
            % Add a small amount to the magnitude to ensure that it's never 0.
            magnitude = norm(tmp_blk_Hist{blk_idx}(:)) + 0.01;
            
            % Divide all of the histogram values by the magnitude to normalize 
            % them.
            tmp_blk_HistN{blk_idx} = tmp_blk_Hist{blk_idx} / magnitude;
            
            % Length of HoG descriptor would increase if we simply
            % concatenate HoG descriptor from member cells here, thus we
            % need to sum gradient component in each direction bin to keep
            % feature size consistency in later HoG contrast calculation.
            for i = 1:hog.numBins
                tmp_blk_HistN_Bin(blk_idx, i) = sum(tmp_blk_HistN{blk_idx}(:,:,i), 'all');
            end
        end
        
        % Calculate local HoG saliency with manner given in formula (7) and
        % (8) given in section 3.2
        % Normalize HoG feature of Target cell 
        THoG_S(1,:) = tmp_Hist(hog.Ty_c, hog.Tx_c, :);  % reshape target cell HoG feature into vector
        mag_THoG = norm(THoG_S) + 0.01;
        T_HoGN = THoG_S(1,:) / mag_THoG;
        
        % calculate HoG contrast between Target cell and local contrast
        % blocks expanded in current sample window
        diff_HoG = zeros(hog.nP, hog.numBins);  % HoG difference vector 
        diff_HoG_scalar = zeros(1, hog.nP);     % HoG difference volume
        for i = 1:hog.nP
            % Difference HoG vector between kernel cell T and i-th local
            % contrast block
            diff_HoG(i, :) = tmp_blk_HistN_Bin(i, :) - T_HoGN;
            diff_HoG_scalar(i) = norm(diff_HoG(i, :),2);    % formula (8)
        end
        % average local DHOGM within 'hog.nP' member cells in local contrast block 
        avg_DHOGM = sum(diff_HoG_scalar(:))/hog.nP; % formula (10)
        HdT = mag_THoG - avg_DHOGM; % formula (9)
        % update HoG saliency map
        HoG_Saliency(y,x) = HdT;
    end
end

% calculate final saliency
out = zeros(imgR, imgC);
for y = 1:imgR
    for x = 1:imgC
        out(y,x) = Sa(y,x) * HoG_Saliency(y,x);
    end
end
out = abs(out);
bw = out_bw(out, 3);

% HoG特征的运算时间较长，对复杂地物背景滤波结果仍有较多虚警
subplot(2,2,1);
imshow(spatial_saliency,[]);
subplot(2,2,2);
imshow(Sa,[]);
subplot(2,2,3);
imshow(HoG_Saliency,[]);
subplot(2,2,4);
imshow(bw,[]);
end

