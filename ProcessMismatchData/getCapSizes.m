function [capSizeDAC, capSizeMAC] = getCapSizes(NOB,process_param, unitCap, capRatio, tolerableGainError, F_per_M)

unitArea = unitCap/F_per_M;
sigmaUnit = process_param/sqrt(unitArea);

LSB = 1/2^(NOB-1);

%% DAC
sum=0;
for i=0:NOB-2
    sum = sum + (2^i);
end


%INL
minAreaINL = (6*process_param*sqrt(sum)/LSB)^2;
DACSizeMultiplier1 = minAreaINL/unitArea;
capSizeDAC1 = ceil(DACSizeMultiplier1)*unitCap
%DNL
sum = sum + 2^(NOB-1);
minAreaINL = (3*process_param*sqrt(sum)/LSB)^2;
DACSizeMultiplier2 = minAreaINL/unitArea;
capSizeDAC2 = ceil(DACSizeMultiplier2)*unitCap

capSizeDAC = max(capSizeDAC1, capSizeDAC2);

% m = 1;
% while (3*sigmaMaxINL)>LSB/2
%     m = m + 1;
%     unitArea = m*unitCap/F_per_M;
%     sigmaUnit = process_param/sqrt(unitArea);
%     sigmaMaxINL = sqrt(sum)*sigmaUnit;
% end
% capSizeDAC = m*unitCap

%% MAC
requiredSigma = tolerableGainError/3 * unitCap;
MACSizeMultiplier = ((process_param/requiredSigma)^2)/unitArea;
capSizeMAC = ceil(MACSizeMultiplier)*unitCap;


end

