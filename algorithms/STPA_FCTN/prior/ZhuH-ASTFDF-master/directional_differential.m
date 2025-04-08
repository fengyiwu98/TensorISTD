function [res] = directional_differential(imin, fvec)
%DIRECTIONAL_DIFFERENTIAL
%

[fvr, fvc] = size(fvec);
[imgR ,imgC, dimension]=size(imin);

res = filter2(fvec, double(imin));
if(fvc>1)  % filter in 1*c
    res(1:imgR,1) = 0;
    res(1:imgR,imgC) = 0;
end
if(fvr>1)  % filter in r*1 
    res(1,1:imgC) = 0;
    res(imgR,1:imgC) = 0;
end
end

