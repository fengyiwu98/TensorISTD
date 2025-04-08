function [x,y,RC,run_time]=RTRCN(tenD,H,lambda2,miu)
%% initialize parameters
N=ndims(tenD);
J=size(tenD);
x=rand(J);
y=zeros(J);
No=zeros(J);
L=ceil(N/2);
l=cell(L,1);

z=cell(L,1);
w=zeros(J);
lambda1=0;
sk=zeros(L,1);
SR=1;
for n=1:L
    l{n}=x;
    z{n}=zeros(J);
    order=[n:N 1:n-1];
    lambda1=lambda1+H/sqrt(SR*max([prod(J(order(1:L))),prod(J(order(L+1:N)))]));
    sk(n)=min([prod(J(order(1:L))),prod(J(order(L+1:N)))]);
end
maxiter=200;
% [n1,n2,n3] = size(tenD);
% lambda= 10/ sqrt(max(n1,n2)*n3);  
[n1,n2,n3,n4] = size(tenD);
% lambda= 1/ sqrt(max(max(max(n1,n2),n3),n4));  
l_cs=zeros(J);
z_cs=zeros(J);
RC=nan(maxiter,2);
epsilon=1e-3;
tol=1e-7;
weight=1./L*ones(1,L);
%% ADMM algorithm
t=cputime;
for i=1:maxiter
    % update l^(n) (auxiliary variables of low-rank part)
    preT=sum(y(:)>0);
    for n=1:L
        order=[n:N 1:n-1];
        m=permute(x-z{n}/miu,order);
        M=reshape(m,prod(J(order(1:L))),[]);
        [M,sk(n)]=shrink_matrix(M,weight(n)/miu,sk(n),epsilon,false);
        m=reshape(M,J(order));
        l{n}=ipermute(m,order);
        l_cs=l_cs+l{n};
        z_cs=z_cs+z{n};
    end
    % update x (low-rank part)
    x=(l_cs+z_cs/miu+(tenD-y-No-w/miu))./(L+1);
    l_cs=zeros(J);
    z_cs=zeros(J);
    % update y (sparse part)
    y=shrink_vector((tenD-x-No-w/miu),lambda1/miu);
    % update z^(n)
    for n=1:L
        z{n}=z{n}+miu*(l{n}-x);
    end
    % update w
    w=w+miu*(x+y+No-tenD);
    %update N
    No=miu*(tenD-x-y+w/miu)/(miu+2*lambda2);
    % evaluate recovery accuracy
    miu=min(miu*1.2,1e10);
   currT = sum(y(:) > 0);
    errList    = norm(tenD(:)-x(:)-y(:)-No(:)) / norm(tenD(:));
%    fprintf('tensor ring iterations = %d   difference=%f\n', i, errList);
disp(['iter ' num2str(i)  ...
                   ', err=' num2str(errList)...
                    ',|T|0 = ' num2str(currT)]); 
                
    if errList < tol %||(preT>0 && currT>0 && preT == currT)
        break;  
    end 

end

end