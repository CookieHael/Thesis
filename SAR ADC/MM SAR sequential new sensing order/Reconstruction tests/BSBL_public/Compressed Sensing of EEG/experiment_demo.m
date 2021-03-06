% Experiment demo for the first EEG experiment in: 
%   Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
%   Compressed Sensing of EEG for Wireless Telemonitoring with Low Energy 
%   Consumption and Inexpensive Hardware, accepted by IEEE Trans. on 
%   Biomedical Engineering, 2012
%
% This demo file shows the ability of BSBL-BO to recover signals which are
% non-sparse in the time domain or in any transformed domain
%
% Author: Zhilin Zhang (zhangzlacademy@gmail.com)
% Date  : Sep 12, 2012



clear; close all;

% load the first channel of the common EEG dataset in the EEGLab
load EEGdata_ch1;
% Set loop variables
%K_array = 1:1:20;
%C_factor_array = 0:5:50;
%C_factor_array(1) = 1;
K_array = 1;
C_factor_array = 15;
% load the randomly generated sparse binary matrix
% (you can use other sensing matrices)
%load Phi;             
                     
% epoch length of the dataset is 384
N = 384;
% each epoch data is compressed to 192 data points (compression 50%)
M = 20;

results_meanMSE = zeros(1,100);
results_meanSSIM = zeros(1,100);
results_meanTime = zeros(1,100);
% results_totalCap = zeros(1,100);
bestPhi = 0;

tStart = tic;
for i = 1:length(K_array)
    for j=1:length(C_factor_array)
        for abc = 1:100
        abc
        % Number of ones per sensing matrix column
        K=K_array(i);
        % Scaling factor of sense cap VS MAC holding cap
        C_factor = C_factor_array(j);
        totalC = K*1 + M*C_factor;
%         results_totalCap(i,j) = totalC;
        
         Phi = generateCorrectedSRBM(M,N,K,C_factor);
%         lowest_mse = 1000;
        %mat = load("Good matrices.mat").remember_this;  %%Capacitor ratio was 35 I think
       % Phi = mat(:,385:768);
%        Phi=load("Phi50ch384MAC.mat").phi;
%         Phi = load("bestPhi.mat").bestPhi;
%         Phi = load("Phi.mat").Phi;
%         Phi = load("sensing192").corrected_Phi;
        
        ADC_count = ceil(9*M/N);
        %fprintf("This needs " + ADC_count + " ADC's or a single one at this overclocking ratio.\n");



        % We use DCT dictionary matrix
        A=zeros(M,N);
        for k=1:M
        A(k,:)=dct(Phi(k,:));
        end

        recon_erp = []; 
        count = 0;
        for ch = 1 : 1       % channel number
            for ep = 1 : 80     % epoch number
                count = count + 1;

                x = EEGdata_ch1(ch,1:N,ep)';
               
                % compress an epoch
                y = Phi * double(x);
 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             use  BSBL-BO to recover
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                blkLen = 13;                  % the block size in the user-defined block partition
                groupStartLoc = 1:blkLen:N;   % user-defined block partition

                % run the algorithm
                tic;
                Result1 = BSBL_BO(A,y,groupStartLoc,0,'prune_gamma',-1, 'max_iters',10);
                t_bsbl(count) = toc;

                % restore recovered data
                recon_EEG(ch,1:N,ep) = (idct(Result1.x))';  

                % MSE
                mse_bsbl(count) = (norm(x - recon_EEG(ch,1:N,ep)','fro')/norm(x,'fro'))^2;
                %rms_bsbl(count) = calculateRMS(x/(max(x)-min(x)), (recon_EEG(ch,1:N,ep)/(max(recon_EEG(ch,1:N,ep))-min(recon_EEG(ch,1:N,ep))))');
                % SSIM
                windowLen = 100;
                [mssim, ssim_map] = ssim_1d( x, recon_EEG(ch,1:N,ep)', windowLen);
                ssim_bsbl(count) = mssim;

                %fprintf('BSBL-BO: time: %4.3f, MSE: %5.4f | SSIM: %5.4f\n',t_bsbl(count),mse_bsbl(count),ssim_bsbl(count));
                
            end
            %fprintf("Average MSE: " + mean(mse_bsbl) + ".\n");
        end
%         filename = "Reconstruction"+K+"_"+C_factor;
%         save(filename,"recon_EEG", "mse_bsbl", "ssim_bsbl", "t_bsbl");
        if (mean(ssim_bsbl)>max(results_meanSSIM))
            bestPhi = Phi;
        end

        results_meanMSE(1,abc) = mean(mse_bsbl);
        results_meanSSIM(1,abc) = mean(ssim_bsbl);
        results_meanTime(1, abc) = mean(t_bsbl);


        % plot the recovery result of an epoch EEG 
%         kepoch = 40;   
%         figure;
%         subplot(411); plot(EEGdata_ch1(1,1:N,kepoch),'linewidth',1); 
%         title('An Original EEG Epoch');
%         subplot(412); plot(recon_EEG(1,:,kepoch),'r','linewidth',1);
%         title_text = ['The Recovered EEG. MSE=',num2str(mse_bsbl(kepoch)), '. SSIM=', num2str(ssim_bsbl(kepoch))];
%         title(title_text);
%         z = max(mse_bsbl);
%         kepoch = find(mse_bsbl==z);
%         subplot(413); plot(EEGdata_ch1(1,1:N,kepoch),'linewidth',1); 
%         title('An Original EEG Epoch');
%         subplot(414); plot(recon_EEG(1,:,kepoch),'r','linewidth',1);
%         title_text = ['The Recovered EEG. MSE=',num2str(mse_bsbl(kepoch)), '. SSIM=', num2str(ssim_bsbl(kepoch))];
%         title(title_text);
%         mn = mean(mse_bsbl);
%           if (mn<lowest_mse)
%               lowest_mse = mn;
%               bestPhi = Phi;
%           end
        end
        %save recovery_result_by_BSBL;
    end
end

totalTime = toc(tStart);
% save("SimResults", "results_meanMSE", "results_meanSSIM", "results_meanTime", "results_totalCap");

