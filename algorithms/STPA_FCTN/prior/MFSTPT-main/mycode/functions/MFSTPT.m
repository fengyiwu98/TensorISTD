function [L,S] = MFSTPT(X, lambda, tenW, opts)
tol =5*1e-4; 
max_iter = 250;
rho = opts.rho;
mu = opts.mu;
max_mu = 1e10;
DEBUG = 1;
N = rankN(X,0.1);%模仿pstnn区秩的方式
epsilon=opts.epsilon;
if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
if isfield(opts, 'N');           N = opts.N;                  end

dim = size(X);
L = zeros(dim);                   
%L = gpuArray(single(L));
S = zeros(dim);
%S = gpuArray(single(S));
Y = zeros(dim);
%Y = gpuArray(single(Y));
weightTen = ones(dim);
%weightTen = gpuArray(single(weightTen));
for iter = 1 : max_iter
    preT = sum(S(:) > 0);
    % update L=B
    R = -S+X-Y/mu;%X=D
    L =Laplace(R, 1/mu, L,epsilon);%
    % update S=T
    T = -L+X-Y/mu;
    %S1 = prox_l1(T1, weightTen*lambda/mu);  
    S = prox_l1(T, weightTen*lambda/mu); %S=tk+1
    weightTen = N./ (abs(S) + 0.01)./tenW;
    dY = L+S-X;
    %dY = -L-S+X;
    err = norm(dY(:))/norm(X(:));
    if DEBUG
        if mod(iter,15) == 0            
            disp(['iter=' num2str(iter) ...
                   ', err=' num2str(err)])... 
                   % ',|T|0 = ' num2str(sum(S(:) > 0))]); 
        end
    end
    currT = sum(S(:) > 0);
    if err < tol || (preT>0 && currT>0 && preT == currT)
        break;
    end 
    Y = Y + dY*mu;
    mu = min(rho*mu,max_mu);    
end

function N= rankN(X, ratioN)
    [n1,~,n3] = size(X);
    D = Unfold(X,n3,1);
    [~, S, ~] = svd(D, 'econ');
    [desS, ~] = sort(diag(S), 'descend');
    ratioVec = desS / desS(1);
    idxArr = find(ratioVec < ratioN);
    if isempty(idxArr)
        N=n1;
    else
        if idxArr(1) > 1
            N = idxArr(1) - 1;
        else 
            N = 1;
        end
    end
    if N>n1
        N=n1;
    end  
