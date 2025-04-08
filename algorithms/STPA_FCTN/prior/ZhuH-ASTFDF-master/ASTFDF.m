function [result_video] = ASTFDF(path, ext, ck, tspan,seq)
%ASTFDF -
% 'Infrared moving point target detection based on an anisotropic 
% spatial-temporal fourth-order diffusion filter' by Zhu H et al in 2018.
% params:
%   'ck' - contrast parameter given in Eq.2 for diffusivity function.
%   'tspan' - temporal evaluation span

opt_res = false; % output threshold segentation result into video
test = false;
vis = true;
fidx = 1;
data=[num2str(seq) '\'];
result_path=['D:\code\tensor_code\shiyan\shiyan\comparision\FCTN_prior_2\prior_result\ASTFDF\sequence',data];
%% set path
addpath(path);
img_dir = fullfile([path '/*.' ext]);

% Create the operators for computing image derivative at every pixel.
hx = [1,0,-1]./2;
hy = hx';
hxx = [1 -2 1];
hyy = hxx';

% init image sequence/video consts
frame_list=dir(img_dir);
n_frm = length(frame_list); % total frame quantity

im0 = imread(frame_list(1).name);
[imgR, imgC, dim] = size(im0);
frame_data = zeros(imgR, imgC, n_frm);
pred = zeros(imgR, imgC, n_frm);
imout = zeros(imgR, imgC, n_frm);
%frame_data(:,:,1) = im0(:,:,1);

% load all image files and compute temporal AFDF
% dx = zeros(imgR, imgC, n_frm);
% dy = zeros(imgR, imgC, n_frm);
for k = 1:(n_frm - 1)
    %frame_data(:,:,k) = rgb2gray(imread(frame_list(k).name));
    frame_data(:,:,k) = im2double(rgb2gray(imread([num2str(k-1) '.' ext])));
    
    % Compute the 1st rank derivative alone x and y direction
    dx = directional_differential(frame_data(:,:,k), hx);
    dy = directional_differential(frame_data(:,:,k), hy);
    
    % Convert the gradient vectors to polar coordinates (angle and magnitude).
    dxs = dx.^2;
    dys = dy.^2;
    dxpdy = dx.*dy;
    magnit = (dys + dxs).^.5;   % L2-norm of gradient of u in discrete domain defined in Eq.20
    
    % compute 2nd rank derivative alone x and y direction
    dxx = directional_differential(frame_data(:,:,k), hxx);
    dyy = directional_differential(frame_data(:,:,k), hyy);
    dxy = directional_differential(frame_data(:,:,k), [1/4 0 -1/4; 0 0 0; -1/4 0 1/4]);
    
    % compute Eta and Epsilon map according to Eq.13~15 in Hajiaboli10's work
    Eta = (dxx.*dxs+2*dxpdy.*dxy+dyy.*dys)./magnit;
    Epsilon = (dxx.*dys-2*dxpdy.*dxy+dyy.*dxs)./magnit;
    
    % compute 2-order temporal gradient
    tgd = zeros(imgR, imgC);
    if(k>1)
        tgd = tgd + frame_data(:,:,k-1) - frame_data(:,:,k);
    end
    if(k<n_frm)
        tgd = tgd + frame_data(:,:,k+1) - frame_data(:,:,k);
    end
    
    % compute temporal differential 'dt' given in Eq.10
    diff = ck^2./(ck^2+magnit.^2);  % diffusivity function defined by Eq.2
    sum_var = diff.^2.*Eta + diff.*Epsilon + tgd;   % Eq.10
    sum_var_dxx = directional_differential(sum_var, hxx);
    sum_var_dyy = directional_differential(sum_var, hyy);
    
    % prediction frame given by Eq.21
    % Eq.21 defines a vector sum operation based on partial derivative map
    % computed from Eq.19, whose prototype can be abstracted as 
    %   'g(x,y,t+k) = g(x,y,t) + dt*[d(dg/dt)/dx^2 d(dg/dt)/dy^2]'
    % the 'dg/dt' is the 'sum_var' afore-computed, and so the laplacian
    % derivatives as update components.

    pred(:,:,k + 1) = frame_data(:,:,k) - tspan.*(sum_var_dxx+sum_var_dyy+tgd);    
    imout(:,:,k + 1) = imGrayNorm(frame_data(:,:,k + 1) - pred(:,:,k + 1), true);
    
    % visualization
    if test
        figure(1)
        subplot(1,3,1);
        imshow(dx,[]);
        subplot(1,3,2);
        imshow(dy,[]);
        subplot(1,3,3);
        imshow(magnit,[]);

        figure(2);
        subplot(1,3,1);
        imshow(imin,[]);
        subplot(1,3,2);
        imshow(Eta,[]);
        subplot(1,3,3);
        imshow(Epsilon,[]);
    end
    
    if vis
        fhnd = figure(fidx);
        set(fhnd,'name',['frame' num2str(k) '_ck-' num2str(ck) '_tspan-' num2str(tspan)]);
        subplot(1,3,1);
        imshow(frame_data(:,:,k),[]);
        subplot(1,3,2);
        imshow(pred(:,:,k + 1),[]);
        subplot(1,3,3);
        imshow(imout(:,:,k + 1),[]);
    end
end

% if opt_res
%     mkdir('results');
%     result_video = img2video(imout, n_frm, [pwd '/results/'],[datestr(now, 30) '_ck-' num2str(ck) '_tspan-' num2str(tspan)],'.avi',25);
% end
if ~exist(result_path, 'dir')
    mkdir(result_path);
end
if opt_res
    img2video(bw, n_frm, [pwd '/results/'],[datestr(now, 30) '_k-' num2str(segK)],'.avi',25);
else
    imwrite(bw(:,:,k),[result_path,num2str(k,'%04d'), '.png']);
        
end
end

