function [RGF] = genRotationGaussianFilter(sig1, sig2, rad, rank, use_softmax)
%GENROTATIONGAUSSIANFILTER
% Generate a gaussian filter of size (2*rank + 1)*(2*rank + 1) with local
% rotation transform mat defined by orientation offset 'rad'.
% This function now demands 2 sigma values for 2D filter matrix.
% Covariance of X and Y under norm-distrib is ALWAYS 0, and same for their
% derivatives.

len = rank * 2 + 1;
RGF = zeros(len);
fac1 = 2*sig1^2;
fac2 = 2*sig2^2;
rot_mat = [cos(rad) sin(rad); -sin(rad) cos(rad)];
%stride = .1;
pos_x(1, 1:len) = -rank:1:rank;
pos_y(1, 1:len) = -rank:1:rank;


for i = 1:len
    for j = 1:len
        trns = rot_mat * [pos_y(i); pos_x(j)];
        RGF(i,j) = exp(-(trns(1)^2/fac2 + trns(2)^2/fac1));
    end
end


if use_softmax
    % softmax origin RGF to compress the variance of sample weights
    soft_RGF = exp(RGF(:,:));
    sum_soft = sum(sum(soft_RGF(:,:)));
    RGF = soft_RGF/sum_soft;
else
    % use simple normalization
    RGF = RGF / sum(sum(RGF));
end
%RGF = RGF/min(min(RGF));
% %surf(pos_x, pos_y, RGF);
end

