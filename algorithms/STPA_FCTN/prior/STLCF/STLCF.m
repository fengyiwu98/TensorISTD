function [] = STLCF(path, ext, tspan_rng, swind_rng, stlc_k,seq)
%STLCF 
% A matlab realisation on 
% Infrared moving point target detection based on an anisotropic spatial-temporal fourth-order diffusion filter'
% by 'Deng L, Zhu H, Tao C, Wei Y' in 'International Journal of Optics.
% 2019 Dec'.
%
% param:
%   'tspan_rng' - time span range for TLC sampling, refers param 'Nt/2' of
%   Eq.2
%   'swind_rng' - spatial span range for SLC sampling, the scale of image
%   window.

opt_res = false; % output threshold segentation result into video
data=[num2str(seq) '\'];
result_path=['D:\code\tensor_code\shiyan\comparision\FCTN_prior_2\prior_result\STLCF\sequence',data];

%% set path
addpath(path);
img_dir = fullfile([path '/*.' ext]);

% init image sequence/video consts
frame_list=dir(img_dir);
n_frm = length(frame_list); % total frame quantity
%n_frm = 30;

im0 = imread(frame_list(1).name);
[imgR, imgC, dim] = size(im0);
frame_data = zeros(imgR, imgC, n_frm);
%frame_data(:,:,1) = im0(:,:,1);
% load all image files
for k = 1:n_frm
    %frame_data(:,:,k) = rgb2gray(imread(frame_list(k).name));
    frame_data(:,:,k) = rgb2gray(imread([num2str(k,'%04d') '.' ext]));
end

% compute STLCF from enumerated image sequence
tlc = zeros(imgR, imgC, n_frm); % temporal local contrast
slc = zeros(imgR, imgC, n_frm); % spatial local contrast
It = zeros(imgR, imgC, n_frm); % normalized TLC
Is = zeros(imgR, imgC, n_frm); % normalized SLC

% init spatial sample span
ssmp_y_prm = zeros(imgR, 3);
ssmp_x_prm = zeros(imgC, 3);
for y = 1:imgR
    % spatial sampling consts alone y-axis
    smp_y_min = max(1, y - swind_rng);
    smp_y_max = min(imgR, y + swind_rng);
    smp_y_span = smp_y_max - smp_y_min + 1;
    ssmp_y_prm(y,:) = [smp_y_min, smp_y_max, smp_y_span];
end
for x = 1:imgC
    smp_x_min = max(1, x - swind_rng);
    smp_x_max = min(imgC, x + swind_rng);
    smp_x_span = smp_x_max - smp_x_min + 1;
    ssmp_x_prm(x,:) = [smp_x_min, smp_x_max, smp_x_span];
end

tsmp_prm = zeros(n_frm, 3);
for k = 1:n_frm
    % temporal sampling consts, start and end frame index, and frame span 
    smp_h = max(1,k-tspan_rng);
    smp_t = min(k+tspan_rng,n_frm);
    n_tsmp = smp_t - smp_h + 1; % sample size in temporal window
    tsmp_prm(k,:) = [smp_h, smp_t, n_tsmp];
    
    frm_max = 0;
    for y = 1:imgR
        for x = 1:imgC
            % TLC within neighbor frames
            tlc(y,x,k) = sum(frame_data(y, x, smp_h:smp_t))/n_tsmp;
            % SLC within neighbor pixels
            % MEMO: you may replace SLC with output saliency of ANY viable
            % single frame detector.
            % spatial sampling consts alone x-axis
            nsmp_pix = ssmp_y_prm(y,3) * ssmp_x_prm(x,3) - 1; % total quantity of sampled neighbor pixel, kernel pix excluded
            mean_itns = (sum(frame_data(ssmp_y_prm(y,1):ssmp_y_prm(y,2), ssmp_x_prm(x,1):ssmp_x_prm(x,2), k),'all') - frame_data(y,x,k))/ nsmp_pix;
            slc(y, x, k) = abs(frame_data(y,x,k) - mean_itns);
            frm_max = max(frm_max, slc(y, x, k));
        end
    end
    
    % normalize SLC with manner of 'picking local max SLC from kernel pixel' given by Eq.6 
    % However, test result shows that normalization below even make spatial
    % saliency worse in ground scene. What a ridiculous outcome.
%     for y = 1:imgR
%         for x = 1:imgC
%             % It seems that SLC normalization refers to global max of
%             % current frame instead of local.
%             %Is(y,x,k) = slc(y, x, k)/max(slc(ssmp_y_prm(y,1):ssmp_y_prm(y,2), ssmp_x_prm(x,1):ssmp_x_prm(x,2), k),[],'all');
%             Is(y,x,k) = slc(y, x, k)/max(slc(:, :, k),[],'all');
%         end
%     end
    Is(:,:,k) = slc(:, :, k)/frm_max;
    %Is(:,:,k) = slc(:, :, k);   % let's try to get rid of spatial normalization
    fprintf('%g/%g\r', k, n_frm);
end

% nomalization on TLC and SLC, and compute STLCF
smp_tmp_max = zeros(imgR, imgC);    % temporal local max of frame k
Ist = zeros(imgR, imgC, n_frm); % STLCF
bw = zeros(imgR, imgC, n_frm); % segmentation result
for k = 1:n_frm
    for y = 1:imgR
        for x = 1:imgC
            smp_tmp_max(y, x) = max(tlc(y,x,tsmp_prm(k,1):tsmp_prm(k,2)),[],'all');  % update max temporal saliency of kernel position
            It(y,x,k) = tlc(y,x,k)/smp_tmp_max(y, x);
            Ist(y,x,k) = It(y,x,k) * Is(y,x,k);
        end
    end
    
    bw(:,:,k) = out_bw(Ist(:,:,k), stlc_k);
    % visualization 
    fhnd = figure(1);
    set(fhnd,'name',['frame' num2str(k)]);
    subplot(1,3,1);
    imshow(frame_data(:,:,k),[]);
    subplot(1,3,2);
    imshow(Ist(:,:,k),[]);
    subplot(1,3,3);
    imshow(bw(:,:,k),[]);


if ~exist(result_path, 'dir')
    mkdir(result_path);
end
if opt_res
    img2video(bw, n_frm, [pwd '/results/'],[datestr(now, 30) '_k-' num2str(segK)],'.avi',25);
else
    %imwrite(bw(:,:,k),[result_path, '/',num2str(k,'%04d'), '.png']);
    imwrite(Ist(:,:,k),[result_path, '/',num2str(k,'%04d'), '.png']);
        
end

end
end

