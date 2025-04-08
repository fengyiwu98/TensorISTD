function [E] = prox_l21(X,beta,mu)
[~,~,n3]=size(X);
for i=1:n3
    L=X(:,:,i);
    A = sqrt(sum(sum(L.^2)));
    E(:,:,i) = (1-beta/(mu*A))*L;
end