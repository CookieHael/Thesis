function Result = STSBL_EM(Phi, y, blkStartLoc, LearnLambda, varargin)
% STSBL-EM: Spatio-temporal Sparse Bayesian Learning Based on the EM Method
% Describe:
%     The algorithm solves the following inverse problem 
%                        Y = Phi * X + V
%     where Y is N by L known matrix, Phi is N by M known matrix, 
%     X is M by L unknown matrix, V is N by L unknown noise matrix. Each
%     column in X has the same block structure (defined by input variable
%     'blkStartLoc'); the locations of nonzero blocks in each column are
%     identical. X may have both intra-column correlation and inter-column
%     correlation, which are exploited to improve the algorithm's
%     performance.
%
%
% ============================== INPUTS ============================== 
%   Phi         : N X M known matrix
%
%   y           : N X L measurement vector 
%
%   blkStartLoc : Start location of each block
%   
%   LearnLambda : (1) If LearnLambda = 1, use the lambda learning rule for generaly noisy 
%                     cases (SNR<=20dB) (thus the input lambda is just as initial value)
%                 (2) If LearnLambda = 2, use the lambda learning rule for high SNR cases (SNR>20dB)
%                 (3) If LearnLambda = 0, do not use the lambda learning rule, but use the input 
%                     lambda value (or the default lambda value).
%                 
%
% [varargin values -- in most cases you can use the default values]
%
%  'LEARNTYPE'    : LEARNTYPE = 0: Ignore intra-column correlation (but
%                   still exploits inter-column correlation)
%                   LEARNTYPE = 1: Exploit intra-block correlation 
%                 [ Default: LEARNTYPE = 1 ]
%
%  'PRUNE_GAMMA'  : threshold to prune out small gamma_i 
%                   (generally, 10^{-3} or 10^{-2})
%
%  'LAMBDA'       : user-input value for lambda
%                  [ Default: LAMBDA=1e-14 when LearnLambda=0; LAMBDA=std(y)*1e-2 in noisy cases ]
%
%  'MAX_ITERS'    : Maximum number of iterations.
%                 [ Default value: MAX_ITERS = 800 ]
%
%  'EPSILON'      : Solution accurancy tolerance parameter 
%                 [ Default value: EPSILON = 1e-6   ]
%
%  'PRINT'        : Display flag. If = 1: show output; If = 0: supress output
%                 [ Default value: PRINT = 0        ]
%
% ==============================  OUTPUTS ============================== 
%   Result : A structured data with:
%      Result.x          : the estimated X
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda
%
%
% ================================ Reference =============================
%   [1] Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Zhouyue Pi, 
%       Bhaskar D. Rao, Spatiotemporal Sparse Bayesian Learning with 
%       Applications to Compressed Sensing of Multichannel Physiological 
%       Signals, IEEE Trans. on Neural Systems and Rehabilitation 
%       Engineering,  vol. 22, no. 6, 2014, DOI: 10.1109/TNSRE.2014.2319334
%
% ============= Author =============
%   Zhilin Zhang (zhilinzhang@ieee.org)
%
% ============= Version =============
%   1.0 (03/13/2014) 
%


% Scaling (will be re-scalling back in the end of the code)
% This step is to make those default values effective
scl = mean(std(y));
if (scl < 0.4) | (scl > 1)
    y = y/scl*0.4;
end

% get problem dimension
[N,M] = size(Phi);
L = size(y,2);

% Default Parameter Values for Any Cases
EPSILON       = 1e-6;       % solution accurancy tolerance  
MAX_ITERS     = 800;        % maximum iterations
PRINT         = 0;          % don't show progress information
LEARNTYPE     = 1;          % adaptively estimate the covariance matrix B

if LearnLambda == 0  
    lambda = 1e-12;   
    PRUNE_GAMMA = 1e-3;  MatrixReg = zeros(L);
elseif LearnLambda == 2
    lambda = 1e-3;    
    PRUNE_GAMMA = 1e-2;  MatrixReg = eye(L) * 0.5;
elseif LearnLambda == 1
    lambda = 1e-3;    
    PRUNE_GAMMA = 1e-2;  MatrixReg = eye(L) * 2;
else
    error(['Unrecognized Value for Input Argument ''LearnLambda''']);
end

 

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'learntype'
                LEARNTYPE = varargin{i+1};
                if LEARNTYPE ~= 1 & LEARNTYPE ~= 0
                    error(['Unrecognized Value for Input Argument ''LEARNTYPE''']);
                end
            case 'prune_gamma'
                PRUNE_GAMMA = varargin{i+1}; 
            case 'lambda'
                lambda = varargin{i+1};    
            case 'epsilon'   
                EPSILON = varargin{i+1}; 
            case 'print'    
                PRINT = varargin{i+1}; 
            case 'max_iters'
                MAX_ITERS = varargin{i+1};  
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end



if PRINT
    fprintf('\n====================================================\n');
    fprintf('           Running STSBL_EM....... \n');
    fprintf('           Information about parameters...\n');
    fprintf('====================================================\n');
    fprintf('PRUNE_GAMMA  : %e\n',PRUNE_GAMMA);
    fprintf('lambda       : %e\n',lambda);
    fprintf('LearnLambda  : %d\n',LearnLambda);    
    fprintf('LearnType    : %d\n',LEARNTYPE);
    fprintf('EPSILON      : %e\n',EPSILON);
    fprintf('MAX_ITERS    : %d\n\n',MAX_ITERS);
end


%% Initialization

y0 = y;
Phi0 = Phi;
blkStartLoc0 = blkStartLoc;
p = length(blkStartLoc);   % block number
for k = 1 : p-1
    blkLenList(k) = blkStartLoc(k+1)-blkStartLoc(k);
end
blkLenList(p) = M - blkStartLoc(end)+1;
maxLen = max(blkLenList);
if sum(blkLenList == maxLen) == p, 
    equalSize = 1;
else
    equalSize = 0;
end

for k = 1 : p
    Sigma0{k} = eye(blkLenList(k));
end

gamma = ones(p,1);
keep_list = [1:p]';
usedNum = length(keep_list);
mu_t = zeros(M,L);
count = 0;

B = eye(L);

%% Iteration
while (1)
    count = count + 1;

    %=========== Prune weights as their hyperparameters go to zero ==============
    if (min(gamma) < PRUNE_GAMMA)
        index = find(gamma > PRUNE_GAMMA);
        usedNum = length(index);
        keep_list = keep_list(index);
        if isempty(keep_list), 
            fprintf('\n====================================================================================\n');
            fprintf('x becomes zero vector. The solution may be incorrect. \n');
            fprintf('Current ''prune_gamma'' = %g, and Current ''EPSILON'' = %g.\n',PRUNE_GAMMA,EPSILON);
            fprintf('Try smaller values of ''prune_gamma'' and ''EPSILON'' or normalize ''y'' to unit norm.\n');
            fprintf('====================================================================================\n\n');
            break; 
        end;
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        
        % prune gamma and associated components in Sigma0 
        gamma = gamma(index);
        temp = Sigma0;
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        
        % construct new Phi
        temp = [];
        for k = 1 : usedNum
            temp = [temp, Phi0(:,blkStartLoc(k):blkStartLoc(k)+blkLenList(k)-1)];
        end
        Phi = temp;
        %clear temp;
    end

    % =================== Spatially Whitening ==============
    y = y0 * sqrtm(inv(B));
    
    %=================== Compute new weights =================
    PhiAPhi = zeros(N);
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        PhiAPhi = PhiAPhi + Phi(:, currentSeg)*Sigma0{i}*Phi(:, currentSeg)';
        currentLoc = currentSeg(end);
    end

    H = Phi' /(PhiAPhi + lambda * eye(N));
    Hy = H * y; 
    HPhi = H * Phi;
    
    mu_x = zeros(size(Phi,2),L);
    Sigma_x = [];
    Cov_x = [];
     
    A = []; invA = []; A0 = zeros(maxLen); r0 = zeros(1); r1 = zeros(1);
   
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        seg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        mu_x(seg,:) = Sigma0{i} * Hy(seg,:);       % solution
        Sigma_x{i} = Sigma0{i} - Sigma0{i} * HPhi(seg,seg) * Sigma0{i};
        %Cov_x{i} = Sigma_x{i} + mu_x(seg,:) * mu_x(seg,:)';
        
        currentLoc = seg(end);
              
        
        %=========== Learn intra-block correlation within columns ===========
        % do not consider correlation structure in each block
        if LEARNTYPE == 0
            A{i} = eye(currentLen);
            invA{i} = eye(currentLen);
 
        % constrain all the blocks have the same correlation structure
        elseif LEARNTYPE == 1
            if equalSize == 0
                if currentLen > 1
                    %temp = Cov_x{i}/gamma(i);
                    temp = Sigma_x{i}/gamma(i) + mu_x(seg,:)*mu_x(seg,:)'/L/gamma(i);
                    
                    r0 = r0 + mean(diag(temp));
                    r1 = r1 + mean(diag(temp,1));
                end
            elseif equalSize == 1
                %A0 = A0 + Cov_x{i}/gamma(i);
                A0 = A0 + Sigma_x{i}/gamma(i) + mu_x(seg,:)*mu_x(seg,:)'/L/gamma(i);
            end

        end % end of learnType

    end % end of usedNum
    
    %=========== Learn the intra-block correlation with regularization ===========
    % If blocks have the same size
    if (equalSize == 1) & (LEARNTYPE == 1)

        % Constrain all the blocks have the same correlation structure
        % (an effective strategy to avoid overfitting)
        b = (mean(diag(A0,1))/mean(diag(A0)));
        if abs(b) >= 0.99, b = 0.99*sign(b); end;
        bs = [];
        for j = 1 : maxLen, bs(j) = (b)^(j-1); end;
        A0 = toeplitz(bs);
 
        for i = 1 : usedNum
             
            A{i} = A0;
            invA{i} = inv(A{i});
        end
    
    % if blocks have different sizes
    elseif (equalSize == 0) & (LEARNTYPE == 1)
        r = r1/r0; if abs(r) >= 0.99, r = 0.99*sign(r); end;

        for i = 1 : usedNum
            currentLen = size(Sigma_x{i},1);

            bs = [];
            for j = 1 : currentLen, bs(j) = r^(j-1); end;
            A{i} = toeplitz(bs);
            invA{i} = inv(A{i});
        end 
    end

    
    % estimate gamma(i) and lambda 
    if LearnLambda == 1         % this learning rule is suitable for general noisy cases (SNR<=20dB)
        gamma_old = gamma;
        lambdaComp = 0; 
        currentLoc = 0;
        for i =  1 : usedNum

            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            lambdaComp = lambdaComp + trace(Phi(:,currentSeg)*Sigma_x{i}*Phi(:,currentSeg)');
            currentLoc = currentSeg(end);

            %gamma(i) = trace(invA{i} * Cov_x{i})/currentLen;
            gamma(i) = trace(invA{i}*Sigma_x{i})/currentLen ...
                + sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/(L*currentLen);
            
            Sigma0{i} = A{i} * gamma(i);
            
        end
        lambda = norm(y - Phi*mu_x,'fro')^2/(N*L) + lambdaComp/N;
        
    elseif LearnLambda == 2   % this learning rule is suitable for high SNR cases (>20dB)
        gamma_old = gamma;
        lambdaComp = 0;  
        currentLoc = 0;
        for i =  1 : usedNum       
            lambdaComp = lambdaComp + trace(Sigma_x{i}*invA{i})/gamma(i);
            
            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            
            gamma(i) = trace(invA{i}*Sigma_x{i})/currentLen ...
                + sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/(L*currentLen);
            Sigma0{i} = A{i} * gamma(i);
            
            currentLoc = currentSeg(end);
        end
        lambda = norm(y - Phi * mu_x,'fro')^2/(N*L) + lambda * (length(mu_x) - lambdaComp)/N; 
        
    else   % only estimate gamma(i)
        gamma_old = gamma;
        currentLoc = 0;
        for i =  1 : usedNum
            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            
            gamma(i) = trace(invA{i}*Sigma_x{i})/currentLen ...
                + sum( sum( (mu_x(currentSeg,:)'*invA{i}).*mu_x(currentSeg,:)' ,2))/(L*currentLen);
            
            Sigma0{i} = A{i} * gamma(i);
            
            currentLoc = currentSeg(end);
        end
    end

    
    % ================= Learn Spatial (inter-column) Covariance Matrix B ============
    mu_old = mu_t;
    mu_t = mu_x * sqrtm(B); 
    B = zeros(L); currentLoc = 0;
    for i = 1 : usedNum
      
        currentLoc = currentLoc + 1;
        currentLen = size(Sigma0{i},1);
        seg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        B = B + mu_t(seg,:)' * invA{i} * mu_t(seg,:)/gamma_old(i);  
        % invA should choose the previous value, but it can use the updated
        % value since the value of invA changes slowly.
        
        currentLoc = seg(end);
    end
    B = B/norm(B,'fro');
    B = B + MatrixReg;
    B = B./norm(B);  
    if size(Phi,2) <= size(Phi,1), B = eye(L); end;

    
    % ================= Check stopping conditions, eyc. ==============
        
        % check convergence
        if (size(mu_t) == size(mu_old))
            dmu = max(max(abs(mu_old - mu_t)));
            if (dmu < EPSILON)  break;  end;
        end;
        if (PRINT)
            disp([' iters: ',num2str(count),...
                ' num coeffs: ',num2str(usedNum), ...
                ' min gamma: ', num2str(min(gamma)),...
                ' gamma change: ',num2str(max(abs(gamma - gamma_old))),...
                ' mu change: ', num2str(dmu)]);
        end;
        if (count >= MAX_ITERS), if PRINT, fprintf('Reach max iterations. Stop\n\n'); end; break;  end;
%     end
    
end;



%% Expand hyperparameyers
gamma_used = sort(keep_list);
gamma_est = zeros(p,1);
gamma_est(keep_list,1) = gamma;


%% reconstruct the original signal
x = zeros(M,L);
currentLoc = 0;
for i = 1 : usedNum

    currentLen = size(Sigma0{i},1);
    currentLoc = currentLoc + 1;
    seg = currentLoc : 1 : currentLoc + currentLen - 1;

    realLocs = blkStartLoc0(keep_list(i)) : blkStartLoc0(keep_list(i))+currentLen-1;

    x( realLocs,: ) = mu_t( seg,: );
    currentLoc = seg(end);
end

if (scl < 0.4) | (scl > 1)
    Result.x = x * scl/0.4;
else
    Result.x = x;
end
%Result.x = Result.x * sqrtm(B);

Result.gamma_used = gamma_used;
Result.gamma_est = gamma_est;
Result.A = A;
Result.count = count;
Result.lambda = lambda;
return;

