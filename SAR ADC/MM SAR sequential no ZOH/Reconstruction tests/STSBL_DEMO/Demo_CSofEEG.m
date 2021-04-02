% A simple demo to show how to use STSBL-EM to recover ambulatory EEG
% signals after compression. The EEG signals are raw signals recorded
% during car-driving. So there are huge noise in the signals. The signals 
% are compressed immediately without any filtering, and later are recovered 
% by STSBL-EM.
%
% Author: Zhilin Zhang (zhilinzhang@ieee.org)
%
% Reference: 
%   Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Zhouyue Pi, Bhaskar D. Rao, 
%   Spatiotemporal Sparse Bayesian Learning with Applications to Compressed 
%   Sensing of Multichannel Physiological Signals, IEEE Trans. on Neural 
%   Systems and Rehabilitation Engineering, vol. 22, no. 6, 2014, 
%   DOI: 10.1109/TNSRE.2014.2319334
%


clear;  clc;

% Load data
load data;     % 125 Hz


% generate sensing matrix Phi (80 by 250)
M = 80;
N = 250;

while(1)
    [Phi,flag] = genPhi(M,N, 2);   
    if flag == 1, break; end;
end

 
% Use DCT dictionary matrix
A=zeros(M,N);
for k=1:M
A(k,:)=dct(Phi(k,:));
end



% The user-defined block partition for ST-SBL
blkLen = 25;
groupStartLoc = 1:blkLen:N;


% variable for storing recovered dataset
X_hat = zeros(size(data));
segNb = floor(size(data,2)/N);
 
% run....
for j = 1 : segNb
    
     
    % ==== Compress signal ======
    y = Phi * data(:,(j-1)*N+1:j*N)'; 
    
    
    % ==== Recover signal =====
    % Note 1: Setting 'prune_gamma' to -1 will result in non-sparse solution
    % (i.e. the solution has no zero-blocks). This may be useful when signals
    % are non-sparse or have non-sparse representation, such as ambulatory
    % EEG signals.
    %
    % Note 2: When signals are sparse, or have sparse representation,
    % setting 'prune_gamma' to a small positive value is suitable (e.g.
    % prune_gamma = 1e-5, the optimal value can be found by tunning). In
    % this case, the solution has zero blocks.
    %
    Result = STSBL_EM(A,y,groupStartLoc,0,'prune_gamma',-1, 'max_iters',20);
    
    signal_hat = idct(Result.x);

    
    % store the recovered epoch
    X_hat(:,(j-1)*N+1:j*N) = (signal_hat)';
    
    % calculate the MSE for current epoch
    mse(j) = (norm(data(:,(j-1)*N+1:j*N) - (signal_hat)','fro')/norm(data(:,(j-1)*N+1:j*N),'fro'))^2;
    fprintf('Segment %d out of %d: MSE = %f\n',j,segNb, mse(j));
     
    
end

% show some typical results
figure;plot(data(1,1:2000),'g');hold on;plot(X_hat(1,1:2000),'r');
       legend('Original Signal at Channel 1','Recovered Signal at Channel 1');
 
figure;plot(data(15,1:2000),'g');hold on;plot(X_hat(15,1:2000),'r');
       legend('Original Signal at Channel 15','Recovered Signal at Channel 15');
       
       
 