function [tenB,tenT,tenN,iter] = TNN(tenD,lambda2)
 max_iter=100;
  mu_max = 1e6;

mu=0.05;
dim=size(tenD);
d1=zeros(dim);
N=d1;
T=d1;
rho=1.2;
lambda3=100;
epsilon         = 1e-7;
sig = zeros(min(dim(1),dim(2)),dim(3));

gamma = 4*1e-2;   
iter = 0;

for iter    = 1 : max_iter
   
    B=prox_tnn(tenD-T-N+d1/mu,1/mu);
    %[ B,sig] = DC1(tenD-T-N+d1/mu,mu,sig,gamma);
    
   T=softthre(tenD-B-N+d1/mu,lambda2./mu); 
    
     N          = (mu*(tenD -B -T) + d1)/(mu+2*lambda3);
    
   d1 = d1 + mu*(tenD -B - T-N);

    
    mu= min(rho*mu, mu_max );  
    
    
    %% Stop Criterion 
    errList    = norm(tenD(:)-B(:)-T(:)-N(:)) / norm(tenD(:));
    %fprintf('ASTTV-NTLA: iterations = %d   difference=%f\n', iter, errList);
    change(iter)=(errList);
   
    if errList < epsilon
        break;  
    end 
end

tenB=B;tenT=T;tenN=0;
end
