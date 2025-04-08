function [B,T,RC]=FCTN(D,W)
%function [B,T,RC]=FCTN(D)

%% initialize parameters
H =0.015; 
tol=1e-3;   
maxiter=200;
lambda1=0;

epsilon=1e-3;
miu = 1e-4; 

SR=1;
Nway=size(D);
N=ndims(D);
L=nchoosek(N,2)/2; 
sk=zeros(L,1);
Wmap=ones(Nway);
RC=nan(maxiter,2);
weight=1./L*ones(1,L); 


%% initialization  

B=rand(Nway);
T=zeros(Nway);
E=zeros(Nway);
U=cell(L,1);
C=cell(L,1);
orderL=[1 2 3 4;3 1 2 4;2 3 1 4]; 

for n=1:L
    U{n}=B;
    C{n}=zeros(Nway);
    order=orderL(n,:);
    lambda1=lambda1+H/sqrt(SR*max([prod(Nway(order(1:2))),prod(Nway(order(3:N)))]));
end


%% ADMM algorithm
t=cputime;
for i=1:maxiter  
    preT=sum(T(:)>0);
    
    % update u^(n) (auxiliary variables of low-rank part)
    u_cs=zeros(Nway);
    z_cs=zeros(Nway);
    
    for n=1:L
        order=orderL(n,:);
        m=permute(B-C{n}/miu,order);
        M=reshape(m,prod(Nway(order(1:2))),[]);
        [M,sk(n)]=shrink_matrix(M,weight(n)/miu,sk(n),epsilon,false);
        m=reshape(M,Nway(order));
        U{n}=ipermute(m,order);
        u_cs=u_cs+U{n};
        z_cs=z_cs+C{n};
    end
    
    % update B (low-rank part)
    B=(u_cs+z_cs/miu+(D-T-E/miu))./(L+1);
    u_cs=zeros(Nway);
    z_cs=zeros(Nway);
    
    % update T (sparse part)
    T=shrink_vector((D-B-E/miu),Wmap.*lambda1/miu); 
    %T=shrink_vector((D-B-E/miu),lambda1/miu); %without prior weight
    
    %update W(weignt part)
    Wmap = min( 1./(abs(T)+0.01)./W, 100);

    
    % update C^(n)
    for n=1:L
        C{n}=C{n}+miu*(U{n}-B);
    end
    
    % update E
    E=E+miu*(B+T-D);
    
    %update N
    %No=miu*(D-B-T+E/miu)/(miu+2*lambda2);
   
    % evaluate recovery accuracy 
    miu = min(miu*1.2,1e10);
    currT = sum(T(:) > 0);
    errList = norm(D(:)-B(:)-T(:)) / norm(D(:));
    %fprintf('tensor ring iterations = %d   difference=%f\n', i, errList);
    disp(['iter ' num2str(i)  ...
                   ', err=' num2str(errList)...
                    ',|T|0 = ' num2str(currT)]); 
                
    if errList < tol %||(preT>0 && currT>0 && preT == currT)
        break;  
    end 

end

end