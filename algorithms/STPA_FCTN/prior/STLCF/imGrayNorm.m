function [res] = imGrayNorm(imin, pos_only)
%IMGRAYNORM 此处显示有关此函数的摘要
%   此处显示详细说明
gmin = min(imin,[],'all');
gmax = max(imin,[],'all');
if ~pos_only
    gdel = gmax - gmin + .000001;
    res = (imin - gmin)./gdel;
    fprintf('input image: min/max = %g/%g\n', gmin, gmax);
else
    res = exp(imin)./exp(gmax);
end
end

