function x = softthre(a, tau)
x = max( abs(a) - tau, 0);
%x = sign(a).* max( abs(a) - tau, 0);
end