function [capSizeDAC, capSizeMAC] = getCapSizes(NOB,process_param, unitCap, capRatio, gainError, F_per_M)

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
capSizeDAC1 = ceil(DACSizeMultiplier1)*unitCap;


%DNL
sum = sum + 2^(NOB-1);
minAreaINL = (2*process_param*sqrt(sum)/LSB)^2;
DACSizeMultiplier2 = minAreaINL/unitArea;
capSizeDAC2 = ceil(DACSizeMultiplier2)*unitCap;

capSizeDAC = max(capSizeDAC1, capSizeDAC2);


%% MAC
sigmaCOverC = gainError/(3-(1+gainError)*(1/(1+capRatio))*(3-3/sqrt(capRatio)));
multiplier = ceil((process_param/sigmaCOverC)^2/unitArea);
capSizeMAC = multiplier*unitCap;


end

