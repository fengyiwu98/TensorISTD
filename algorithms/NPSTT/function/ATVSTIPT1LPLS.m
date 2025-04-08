function [tenB, tenT] = ATVSTIPT1LPLS(tenD)
[n1,n2,n3]=size(tenD);
%% produce dicitonary
W = fft(tenD,[],3);
for i=1:n3
    [u,s,v] = svd(W(:,:,i),'econ');
    U(:,:,i)=u;
    S(:,:,i)=s;
    V(:,:,i)=v;
end
U=ifft(U,[],3);S=ifft(S,[],3);V=ifft(V,[],3);
D=tprod(U,S);

%%  Initialization
 [~,n4,~]=size(D);
 X=zeros(n4,n2,n3);
 d1=X;
 A=d1;
d2 = zeros(n1,n2,n3);
T=d2;
lambda=1/(sqrt(n3*max(n1,n2)));
mu=1e-2;
mu_max=1e+7;
tol=1e-7;
ro=1.5;
max_iter=100;
change=zeros(1,max_iter);
DT=tran(D);
D_in=t_inverse(D);
normD=norm(tenD(:));

for iter = 1 : max_iter
  
%% Update A
    A_pre=A;
   
   A=tprod(D_in,X+d1/mu+tprod(DT,tenD-T+d2/mu));

%% Update  X
    X_pre=X;
    X=prox_tnn(A-d1/mu,1/mu);
  % X=prox_tensorLaplace(A-d1/mu,1/mu,X_pre,1e-6);
%% Update T
    T_pre=T;
    M=tprod(D,A);
    T = softthre(tenD-M+d2/mu,lambda/mu);
    
%% check convergence
%% check convergence
%     leq1 = X-A;
%     leq2 = tenD-tprod(D,A)-T;
%     leqm1 = max(abs(leq1(:)));
%     leqm2 = max(abs(leq2(:)));
%     
%     difJ = max(abs(A(:)-A_pre(:)));
%     difE = max(abs(T(:)-T_pre(:)));
%     difZ = max(abs(X(:)-X_pre(:)));
%     err = max([leqm1,leqm2,difJ,difZ,difE]);
%     if ( (iter==1 || mod(iter,20)==0))
%         sparsity=length(find(T~=0));
%         fprintf('iter = %d, err = %.8f, beta=%.2f, sparsity=%d\n'...
%             , iter,err,mu,sparsity);
%     end
%     if err < tol
%         break;
%     end
 
errList    = norm(tenD(:)-M(:)-T(:)) / normD;
   fprintf('ASTTV-NTLA: iterations = %d   difference=%f\n', iter, errList);
    change(iter)=(errList);
    if errList < tol
        break;  
    end 

    
%% Update  d1 d2
    d1 = d1 + mu*( X-A);
    d2 = d2 + mu*(tenD-M-T);
%% Update mu
    mu = min(ro*mu,mu_max);
end
%% Output
tenB=M;tenT=T;