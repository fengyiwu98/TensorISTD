function [L,S] = TRPCA(X, lambda, opts)

    tol = 1e-2; %小于这个数就停止
    max_iter = 500;
    rho = 1.05; 
    mu = 2*1e-3;
    max_mu = 1e10;
    DEBUG = 0;
    N = rankN(X,0.1); %X的秩


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
    S = zeros(dim);
    Y = zeros(dim);

    for iter = 1 : max_iter

        preT = sum(S(:) > 0);
        
        % update L 背景
        R = -S+X-Y/mu;
        L = prox_tnn(R,1/mu);

        % update S 目标
        T = -L+X-Y/mu;
        S = prox_l1(T, lambda/mu);
        
        % updata Y mu 其他参数
        dY = L+S-X;
        err = norm(dY(:))/norm(X(:));

        if DEBUG
            if iter == 1 || mod(iter, 1) == 0
                disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', err=' num2str(err)...
                    ',|T|0 = ' num2str(sum(S(:) > 0))]);
            end
        end

        currT = sum(S(:) > 0);

        if err < tol || (preT>0 && currT>0 && preT == currT)
            break;
        end
        
        Y = Y + dY*mu;
        mu = min(rho*mu,max_mu);
        
    end
end

function N = rankN(X, ratioN)
    [~,~,n3] = size(X);
    D = Unfold(X,n3,1); % �? 
    [~, S, ~] = svd(D, 'econ');
    [desS, ~] = sort(diag(S), 'descend');
    ratioVec = desS / desS(1);
    idxArr = find(ratioVec < ratioN);
    if idxArr(1) > 1
        N = idxArr(1) - 1;
    else
        N = 1;
    end
end
    
function [X] = Unfold( X, dim, i )
     X = reshape(shiftdim(X,i-1), dim(i), []);
end

function x = prox_l1(b,lambda)

x = max(0,b-lambda)+min(0,b+lambda);
x = max(x,0);

end

function [X] = prox_pstnn(Y,N,mu)

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);
tau = 1/mu;


% first frontal slice
[U,S,V] = svd(Y(:,:,1),'econ');
diagS = diag(S);
[desS, sIdx] = sort(diagS, 'descend');
%.[Y1, Y2, Y3, ...] = deal(X1, X2, X3, ...)相当�?Y1 = X1; Y2 = X2; Y3 = X3; ...
[desU, desV] = deal(U(:, sIdx), V(:, sIdx));
[U1, diagS1, V1] = deal(desU(:, 1:N), desS(1:N), desV(:, 1:N));
[U2, diagS2, V2] = deal(desU(:, N+1:end), desS(N+1:end), desV(:, N+1:end));    
threshS2 = max(diagS2-tau, 0);    
X(:,:,1) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';


% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    [U,S,V] = svd(Y(:,:,i),'econ');
    diagS = diag(S);
    [desS, sIdx] = sort(diagS, 'descend');
    [desU, desV] = deal(U(:, sIdx), V(:, sIdx));
    [U1, diagS1, V1] = deal(desU(:, 1:N), desS(1:N), desV(:, 1:N));
    [U2, diagS2, V2] = deal(desU(:, N+1:end), desS(N+1:end), desV(:, N+1:end));    
    threshS2 = max(diagS2-tau, 0);    
    X(:,:,i) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';
    X(:,:,n3+2-i) = conj(X(:,:,i));
end
  
% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    diagS = diag(S);
    [desS, sIdx] = sort(diagS, 'descend');
    [desU, desV] = deal(U(:, sIdx), V(:, sIdx));
    [U1, diagS1, V1] = deal(desU(:, 1:N), desS(1:N), desV(:, 1:N));
    [U2, diagS2, V2] = deal(desU(:, N+1:end), desS(N+1:end), desV(:, N+1:end));    
    threshS2 = max(diagS2-tau, 0);    
    X(:,:,i) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';
end

X = ifft(X,[],3);
end