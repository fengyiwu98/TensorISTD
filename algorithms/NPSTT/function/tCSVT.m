function [L,S] = tCSVT(X,lambda)

    tol = 1e-2; 
    max_iter = 500;  % �?��迭代次数
    rho = 1.05; 
    mu = 1e-3;
    threshold = 3;
    max_mu = 1e10;
    DEBUG = 0;
    
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
    
    % 初始�?
    dim = size(X);  % 时空张量尺寸
    L = zeros(dim);  % 低秩分量B
    S = zeros(dim);  % �?��分量T
    Y = zeros(dim);  % 乘子Z
    
    N = ones(1,1+round(dim(1,3)/2));  % [I_3 / 2] + 1 个元�?
    
    for iter = 1 : max_iter

        preT = sum(S(:) > 0);  % 目标中非0元素的个�?
        
        % update L 背景B
        R =  -S+X-Y/mu;
        L = t_CSVT(R,N,mu);  % TCNN算子作用于R以更新L（N为分段奇异�?个数，tau = 1/mu�?
        % 计算分段奇异值个数（�?��迭代中每个切片的分段位置都不同）
        for index = 1:1+round(dim(1,3)/2)
             [~,S,~] = svd(L(:,:,index),'econ');
             diagS = diag(S);
             thresh = max(diagS-threshold, 0);
             % 分段奇异值个数N(1,index)设为L中大于threshold的奇异�?个数
             N(1,index) = size(find(thresh>0),1);
        end
        
        % update S 目标T
        T = -L+X-Y/mu;
        S = prox_l1(T, lambda/mu);
        
        % updata Y（乘子Z�?mu 其他参数
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

    

function x = prox_l1(b,lambda)

x = max(0,b-lambda)+min(0,b+lambda);
x = max(x,0);

end

function [X] = t_CSVT(Y,N,mu)
% Tensor Capped Singular Value Thresholding(t-CSVT)
% 张量上限奇异值阈值算�?t-CSVT)
% N(1,i)为分段奇异�?个数，tau = 1/mu

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);
tau = 1/mu;

% first frontal slice 第一个正面切�?
% 降序SVD
[U,S,V] = svd(Y(:,:,1),'econ');
diagS = diag(S);
[desS, sIdx] = sort(diagS, 'descend');
%.[Y1, Y2, Y3, ...] = deal(X1, X2, X3, ...)相当�?Y1 = X1; Y2 = X2; Y3 = X3; ...
[desU, desV] = deal(U(:, sIdx), V(:, sIdx));
% 以N为界限，将降序的USV切成前一半和后一�?
[U1, diagS1, V1] = deal(desU(:, 1:N(1,1)), desS(1:N(1,1)), desV(:, 1:N(1,1)));
[U2, diagS2, V2] = deal(desU(:, N(1,1)+1:end), desS(N(1,1)+1:end), desV(:, N(1,1)+1:end));    
threshS2 = max(diagS2-tau, 0);
% 下式等价于[U1,U2]*diag([diagS1,threshS2])*[V1,V2]'
X(:,:,1) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';


% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    [U,S,V] = svd(Y(:,:,i),'econ');
    diagS = diag(S);
    [desS, sIdx] = sort(diagS, 'descend');
    [desU, desV] = deal(U(:, sIdx), V(:, sIdx));
    [U1, diagS1, V1] = deal(desU(:, 1:N(1,i)), desS(1:N(1,i)), desV(:, 1:N(1,i)));
    [U2, diagS2, V2] = deal(desU(:, N(1,i)+1:end), desS(N(1,i)+1:end), desV(:, N(1,i)+1:end));    
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
    [U1, diagS1, V1] = deal(desU(:, 1:N(1,2)), desS(1:N(1,2)), desV(:, 1:N(1,2)));
    [U2, diagS2, V2] = deal(desU(:, N(1,2)+1:end), desS(N(1,2)+1:end), desV(:, N(1,2)+1:end));    
    threshS2 = max(diagS2-tau, 0);    
    X(:,:,i) = U1*diag(diagS1)*V1' + U2*diag(threshS2)*V2';
end

X = ifft(X,[],3);
end
