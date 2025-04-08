function [imout] = AFDF(imin, ck, tspan)
%AFDF - An Anisotropic Fourth-Order Diffusion Filter for Image Noise Removal
% Matlab realisation on image enhance algorithm raised by 'Hajiaboli' in
% his 'An Anisotropic Fourth-Order Diffusion Filter for Image Noise Removal
% codes'.
% This algorithm is also basis of 'Infrared moving point target detection
% based on an anisotropic spatial-temporal fourth-order diffusion filter'
% by Zhu H et al in 2018.
% params:
%   'ck' - contrast parameter given in Eq.2 for diffusivity function.
%   'tspan' - temporal evaluation span

% convert input frame into double format
test = false;
if(ischar(imin))
    imin = imread(imin);
end
[imgR ,imgC, dimension]=size(imin);
if(dimension>2) 
    imin = rgb2gray(imin);
    dimension = 1;
end
imin = im2double(imin);

% Create the operators for computing image derivative at every pixel.
hx = [1,0,-1]./2;
hy = hx';
hxx = [1 -2 1];
hyy = hxx';

% Compute the 1st rank derivative alone x and y direction
dx = directional_differential(imin, hx);
dy = directional_differential(imin, hy);

% Convert the gradient vectors to polar coordinates (angle and magnitude).
angles = atan2(dy, dx);
dxs = dx.^2;
dys = dy.^2;
dxpdy = dx.*dy;
magnit = (dys + dxs).^.5;   % L2-norm of gradient of u in discrete domain defined in Eq.20

% compute 2nd rank derivative alone x and y direction
dxx = directional_differential(imin, hxx);
dyy = directional_differential(imin, hyy);
dxy = directional_differential(imin, [1/4 0 -1/4; 0 0 0; -1/4 0 1/4]);

% compute Eta and Epsilon map according to Eq.13~15 in Hajiaboli10's work
Eta = (dxx.*dxs+2*dxpdy.*dxy+dyy.*dys)./magnit;
Epsilon = (dxx.*dys-2*dxpdy.*dxy+dyy.*dxs)./magnit;

% compute temporal differential 'dt' given in Eq.19
diff = ck^2./(ck^2+magnit.^2);  % diffusivity function defined by Eq.2
sum_var = diff.^2.*Eta + diff.*Epsilon;
sum_var_dxx = -directional_differential(sum_var, hxx);
sum_var_dyy = -directional_differential(sum_var, hyy);
%dut = - (sum_var_dxx.^2 + sum_var_dyy.^2).^.5;  % magnitude of temporal partial derivative in Eq.19

% prediction frame given by Eq.21
% Eq.21 defines a vector sum operation based on partial derivative map
% computed from Eq.19, whose prototype can be abstracted as 
%   'g(x,y,t+k) = g(x,y,t) + dt*[d(dg/dt)/dx^2 d(dg/dt)/dy^2]'
% the 'dg/dt' is the 'sum_var' afore-computed, and so the laplacian
% derivatives as update components.

imout = imin - tspan.*(sum_var_dxx+sum_var_dyy);    

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

figure(1)
subplot(1,3,1);
imshow(imin,[]);
subplot(1,3,2);
imshow(imout,[]);
subplot(1,3,3);
imshow(imin - imout);

end

