function [smp_ofs] = genSmpOfs(k)
%GENSMPOFS 
% Generate a 2-page-k-rank square contains image plane sample offset.
% Page 1 and 2 are for row and col offset respectively.
% 

len = 2*k + 1;
smp_ofs = zeros(len, len, 2);

del = -k;
for y = 1:len
    smp_ofs(y,:,1) = del;
    smp_ofs(:,y,2) = del;
    del = del + 1;
end
end

